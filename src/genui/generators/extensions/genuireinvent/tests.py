# genuireinvent/tests.py
import os
import shutil
import tempfile

from django.apps import apps
from django.conf import settings
from django.contrib.auth import get_user_model
from django.test import override_settings
from django.urls import reverse

from rest_framework import status
from rest_framework.test import APITestCase

from genui.qsar.tests import QSARModelInit
from genui.models.models import Algorithm, AlgorithmMode, ModelFileFormat

from . import models
from .tasks import runReinventStagedLearning


TEST_EPOCHS = 2


# ---------------------------------------------------------------------
# Hard dependency checks (NO fakes, NO dummy processes)
# ---------------------------------------------------------------------
def _get_reinvent_bin() -> str | None:
    return (
        getattr(settings, "REINVENT_BIN", None)
        or os.environ.get("REINVENT_BIN")
        or shutil.which("reinvent")
    )


def _assert_reinvent_runtime_available():
    # 1) REINVENT python package (needed for corpus preprocess)
    try:
        from reinvent.datapipeline import preprocess  # noqa: F401
    except Exception as e:
        raise AssertionError(
            "Missing Python dependency: reinvent.datapipeline.preprocess.\n"
            "Install the REINVENT python package (or ensure it is importable in the test environment)."
        ) from e

    # 2) Prior file
    try:
        prior_path = models._resolve_reinvent_prior_path()
    except FileNotFoundError as e:
        raise AssertionError(
            f"{e}\n\n"
            "Fix by either:\n"
            "  - setting REINVENT_PRIOR (env), or\n"
            "  - setting settings.REINVENT_PRIOR_PATH (or GENUI_SETTINGS['REINVENT_PRIOR_PATH']), or\n"
            "  - placing the prior under GENUI_SETTINGS['FILES_DIR']/checkpoints/prior/reinvent.prior\n"
        ) from e

    if not os.path.isfile(prior_path):
        raise AssertionError(f"Resolved REINVENT prior path does not exist: {prior_path}")

    # 3) REINVENT CLI binary
    reinvent_bin = _get_reinvent_bin()
    if not reinvent_bin:
        raise AssertionError(
            "REINVENT CLI binary not found.\n"
            "Fix by either:\n"
            "  - setting settings.REINVENT_BIN, or\n"
            "  - setting REINVENT_BIN env var, or\n"
            "  - ensuring `reinvent` is on PATH."
        )
    if not os.path.isfile(reinvent_bin) and shutil.which(reinvent_bin) is None:
        raise AssertionError(f"REINVENT bin was set but not found on disk/PATH: {reinvent_bin}")

    # If it's a file, ensure executable
    if os.path.isfile(reinvent_bin) and not os.access(reinvent_bin, os.X_OK):
        raise AssertionError(f"REINVENT bin exists but is not executable: {reinvent_bin}")


# ---------------------------------------------------------------------
# Shared setup + helpers (schema-introspective, post-migrations friendly)
# ---------------------------------------------------------------------
class SetUpReinventMixIn(QSARModelInit):
    """
    Setup:
      - creates admin
      - ensures Algorithm + Mode metadata exist
      - uses a temp MEDIA_ROOT so tests don't pollute repo
      - asserts real REINVENT runtime exists (no fakes)
    """

    @classmethod
    def setUpTestData(cls):
        User = get_user_model()
        cls._admin = User.objects.create_superuser(
            username="fadeevartem",
            email="fadeev19190@gmail.com",
            password="1234",
        )

        cls.mode_generator = AlgorithmMode.objects.get_or_create(name="generator")[0]
        cls.alg_reinvent = Algorithm.objects.get_or_create(name="ReinventNet")[0]
        cls.alg_reinvent.validModes.add(cls.mode_generator)
        cls.alg_reinvent.corePackage = "genui.generators.extensions.genuireinvent.genuimodels"
        cls.alg_reinvent.save(update_fields=["corePackage"])

        fmt, _ = ModelFileFormat.objects.get_or_create(
            fileExtension=".pkg",
            defaults={"description": "State of a neural network built with pytorch."},
        )
        if fmt not in cls.alg_reinvent.fileFormats.all():
            cls.alg_reinvent.fileFormats.add(fmt)

    def setUp(self):
        super().setUp()
        self.client.force_login(self._admin)

        # ensure project ownership for queryset visibility
        if getattr(self.project, "owner_id", None) != self._admin.id:
            self.project.owner = self._admin
            self.project.save(update_fields=["owner"])

        # temp MEDIA_ROOT to avoid repo pollution
        self._tmp_media = tempfile.mkdtemp(prefix="reinvent_test_media_")
        self.addCleanup(lambda: shutil.rmtree(self._tmp_media, ignore_errors=True))

        old_media_root = getattr(settings, "MEDIA_ROOT", None)
        self.addCleanup(lambda: setattr(settings, "MEDIA_ROOT", old_media_root))
        settings.MEDIA_ROOT = self._tmp_media

        # hard-check real runtime (no fakes)
        _assert_reinvent_runtime_available()

    # ----------------- creation helpers (introspection) -----------------
    def _create_with_model_fields(self, model_cls, **kwargs):
        field_names = {f.name for f in model_cls._meta.get_fields()}
        filtered = {k: v for k, v in kwargs.items() if k in field_names}
        return model_cls.objects.create(**filtered)

    def _ensure_dataset_links(self, model_cls, kwargs: dict) -> dict:
        fields = {f.name: f for f in model_cls._meta.get_fields()}
        if "project" in fields and "project" not in kwargs:
            kwargs["project"] = self.project
        if "molecules" in fields and "molecules" not in kwargs:
            kwargs["molecules"] = self.molset
        if "molset" in fields and "molset" not in kwargs:
            kwargs["molset"] = self.molset
        return kwargs

    def _get_builder_model(self):
        for app_label in ("models", "genui_models", "genui"):
            for cls_name in ("Builder", "ModelBuilder"):
                try:
                    return apps.get_model(app_label, cls_name)
                except Exception:
                    continue
        raise RuntimeError("Could not locate Builder model (tried Builder/ModelBuilder).")

    def _create_builder_for(self, *, model_class_name: str):
        Builder = self._get_builder_model()
        kwargs = {}

        for f in Builder._meta.fields:
            if getattr(f, "primary_key", False) or getattr(f, "auto_created", False):
                continue
            if getattr(f, "auto_now", False) or getattr(f, "auto_now_add", False):
                continue
            if getattr(f, "has_default", lambda: False)() and f.has_default():
                continue
            if getattr(f, "null", False) or getattr(f, "blank", False):
                continue

            name = f.name.lower()

            # FKs
            if getattr(f, "many_to_one", False) and getattr(f, "remote_field", None):
                rel = f.remote_field.model
                rel_name = getattr(rel, "__name__", "")
                if rel_name == "Project":
                    kwargs[f.name] = self.project
                    continue
                if rel_name in ("User", get_user_model().__name__):
                    kwargs[f.name] = self._admin
                    continue
                if rel_name == "Algorithm":
                    kwargs[f.name] = self.alg_reinvent
                    continue
                if rel_name == "AlgorithmMode":
                    kwargs[f.name] = self.mode_generator
                    continue
                continue

            # choices
            if getattr(f, "choices", None):
                kwargs[f.name] = f.choices[0][0]
                continue

            internal = f.get_internal_type()
            if internal in ("CharField", "TextField"):
                if "class" in name and "model" in name:
                    kwargs[f.name] = model_class_name
                elif "name" in name:
                    kwargs[f.name] = f"builder:{model_class_name}"
                elif "status" in name or "state" in name:
                    kwargs[f.name] = "created"
                else:
                    kwargs[f.name] = "test"
            elif internal in ("IntegerField", "BigIntegerField", "PositiveIntegerField", "SmallIntegerField"):
                kwargs[f.name] = 0
            elif internal in ("FloatField", "DecimalField"):
                kwargs[f.name] = 0.0
            elif internal == "BooleanField":
                kwargs[f.name] = False
            elif internal == "JSONField":
                kwargs[f.name] = {}
            else:
                kwargs[f.name] = "test"

        return Builder.objects.create(**kwargs)

    def _create_model_like(self, model_cls, **kwargs):
        fields = {f.name: f for f in model_cls._meta.get_fields()}

        if "project" in fields and "project" not in kwargs:
            kwargs["project"] = self.project
        if "algorithm" in fields and "algorithm" not in kwargs:
            kwargs["algorithm"] = self.alg_reinvent
        if "mode" in fields and "mode" not in kwargs:
            kwargs["mode"] = self.mode_generator
        if "builder" in fields and "builder" not in kwargs:
            kwargs["builder"] = self._create_builder_for(model_class_name=model_cls.__name__)

        return self._create_with_model_fields(model_cls, **kwargs)

    def _create_dataset_like(self, model_cls, **kwargs):
        kwargs = self._ensure_dataset_links(model_cls, kwargs)
        return self._create_with_model_fields(model_cls, **kwargs)

    def _create_strategy_like(self, model_cls, *, model_instance, **kwargs):
        fields = {f.name: f for f in model_cls._meta.get_fields()}

        if "modelInstance" in fields and "modelInstance" not in kwargs:
            kwargs["modelInstance"] = model_instance
        if "algorithm" in fields and "algorithm" not in kwargs:
            kwargs["algorithm"] = self.alg_reinvent
        if "mode" in fields and "mode" not in kwargs:
            kwargs["mode"] = self.mode_generator
        if "epochs" in fields and "epochs" not in kwargs:
            kwargs["epochs"] = 1

        return self._create_with_model_fields(model_cls, **kwargs)

    # ----------------- API helpers -----------------
    def _create_reinvent_net_via_api(self, *, build: bool):
        url = reverse("reinvent-net-list")
        payload = {
            "name": "Test Reinvent Network",
            "description": "test",
            "project": self.project.id,
            "build": bool(build),
            "trainingStrategy": {
                "algorithm": Algorithm.objects.get(name="ReinventNet").id,
                "mode": AlgorithmMode.objects.get(name="generator").id,
                "epochs": TEST_EPOCHS,
                "batch_size": 16,
                "sample_batch_size": 100,
                "save_every_n_epochs": 1,
            },
            "validationStrategy": {"validSetSize": 5, "split_method": "random", "valid_fraction": 0.2},
            "molset": self.molset.id,
        }
        resp = self.client.post(url, data=payload, format="json")
        self.assertEqual(resp.status_code, status.HTTP_201_CREATED, msg=resp.data)
        return models.ReinventNet.objects.get(pk=resp.data["id"]), resp.data


# ---------------------------------------------------------------------
# Tests: Transfer learning integration (prepareData + TL subprocess)
# ---------------------------------------------------------------------
@override_settings(
    ROOT_URLCONF="genui.urls",
    CELERY_TASK_ALWAYS_EAGER=True,
    CELERY_TASK_EAGER_PROPAGATES=True,
)
class ReinventTransferLearningIntegrationTests(SetUpReinventMixIn, APITestCase):
    def test_prepare_corpus_endpoint_runs_datapipeline_and_writes_files(self):
        net, _ = self._create_reinvent_net_via_api(build=False)

        url = reverse("reinvent-net-prepare-corpus", kwargs={"pk": net.id})
        resp = self.client.post(url, data={}, format="json")
        self.assertEqual(resp.status_code, status.HTTP_201_CREATED, msg=resp.data)

        # Returned paths must exist
        for key in ("train_file", "valid_file", "preview_file", "full_file"):
            p = resp.data.get(key)
            self.assertTrue(p and os.path.isfile(p), f"{key} missing or not a file: {p}")

        # And counts should be >= 0 (valid may be 0 for very small corpora)
        self.assertGreaterEqual(int(resp.data.get("prepared_train", 0)), 0)
        self.assertGreaterEqual(int(resp.data.get("prepared_valid", 0)), 0)

    def test_transfer_learning_direct_call_runs_and_writes_checkpoint_and_log(self):
        net, _ = self._create_reinvent_net_via_api(build=False)

        # 1) corpus
        train_mf, valid_mf = net.prepareData()
        self.assertTrue(os.path.isfile(train_mf.path))
        self.assertTrue(os.path.isfile(valid_mf.path))
        self.assertTrue(os.path.isfile(net.corpusFullFile.path))
        self.assertTrue(os.path.isfile(net.corpusPreviewFile.path))

        # 2) TL
        ckpt_path = net.run_transfer_learning(device="cpu")
        self.assertEqual(ckpt_path, net.checkpointFile.path)
        self.assertTrue(os.path.isfile(ckpt_path), f"Expected checkpoint at {ckpt_path}")
        self.assertGreater(os.path.getsize(ckpt_path), 0, "Checkpoint file is empty")

        # 3) TOML + log should exist
        self.assertTrue(os.path.isfile(net.tlTomlFile.path))
        self.assertTrue(os.path.isfile(net.trainLogFile.path))
        log_txt = open(net.trainLogFile.path, "r", encoding="utf-8").read()
        self.assertIn("[CMD]", log_txt)

        # 4) best_epoch parsing is optional (depends on REINVENT output)
        ts = net.trainingStrategy
        if getattr(ts, "best_epoch", None) is not None:
            self.assertIsInstance(ts.best_epoch, int)
        if getattr(ts, "best_valid_loss", None) is not None:
            self.assertIsInstance(ts.best_valid_loss, float)

    def test_transfer_learning_via_build_task_executes_full_pipeline(self):
        """
        This exercises the actual build pipeline:
          - ReinventNetViewSet create(build=True) enqueues BuildReinventModel
          - builder.getX() calls prepareData()
          - algorithm.fit() triggers CLI TL
        """
        net, data = self._create_reinvent_net_via_api(build=True)

        # If the viewset returns task_id, eager mode executes immediately anyway.
        # We assert the artifacts exist on disk after the build.
        net.refresh_from_db()

        # Builder pipeline should have created corpus + checkpoint + log
        self.assertTrue(os.path.isfile(net.corpusFullFile.path), "Corpus was not created by build pipeline")
        self.assertTrue(os.path.isfile(net.checkpointFile.path), "Checkpoint was not created by build pipeline")
        self.assertTrue(os.path.isfile(net.trainLogFile.path), "Training log was not created by build pipeline")
        self.assertGreater(os.path.getsize(net.checkpointFile.path), 0, "Checkpoint file is empty")


# ---------------------------------------------------------------------
# Tests: Staged learning integration (TOML build + RL subprocess via Celery)
# ---------------------------------------------------------------------
@override_settings(
    ROOT_URLCONF="genui.urls",
    CELERY_TASK_ALWAYS_EAGER=True,
    CELERY_TASK_EAGER_PROPAGATES=True,
)
class ReinventStagedLearningIntegrationTests(SetUpReinventMixIn, APITestCase):
    def _mk_env_agent_gen(self, net: models.ReinventNet, *, add_diversity: bool = False):
        """Create Environment → AgentTraining → Agent → Reinvent (run).

        Returns (env, agent, gen).
        """
        df = None
        if add_diversity:
            df = models.ReinventDiversityFilter.objects.create(
                type="ScaffoldSimilarity",
                bucket_size=10,
                minscore=0.4,
                minsimilarity=0.4,
                penalty_multiplier=0.5,
            )

        env = self._create_dataset_like(
            models.ReinventEnvironment,
            name="Test RL Environment",
            prior_net=net,
            agent_net=net,
            diversity_filter=df,
            aggregation_type="geometric_mean",
        )

        train_cfg = self._create_strategy_like(
            models.ReinventAgentTraining,
            model_instance=net,
            batch_size=16,
            unique_sequences=True,
            randomize_smiles=True,
            tb_isim=False,
            use_checkpoint=False,
            purge_memories=False,
            summary_csv_prefix="reinvent",
            learning_type="dap",
            sigma=64.0,
            rate=0.0005,
        )

        agent = self._create_model_like(
            models.ReinventAgent,
            name="Test Reinvent Agent",
            description="agent",
            environment=env,
            training=train_cfg,
            validation=None,
            output_model=None,
            tb_logdir=os.path.join(settings.MEDIA_ROOT, "tb_rl"),
            json_out_config="_staged_learning.json",
        )

        gen = self._create_model_like(
            models.Reinvent,
            name="Test Reinvent Run",
            description="run",
            environment=env,
            agent=agent,
        )

        return env, agent, gen

    def _create_property_scorer(self, *, name="QED_scorer", property_name="Qed",
                                weight=1.0, project=None):
        return models.PropertyScorer.objects.create(
            name=name,
            property_name=property_name,
            weight=weight,
            project=project or self.project,
        )

    def _create_stage_with_scorers(self, gen, *, order=0, max_steps=5,
                                    property_scorers=None, aggregation_type="geometric_mean"):
        """Create a ReinventStage and assign PropertyScorers via M2M."""
        stage = models.ReinventStage.objects.create(
            generator=gen,
            order=order,
            termination_type="simple",
            max_score=1.0,
            min_steps=1,
            max_steps=max_steps,
            scoring_source="inline",
            aggregation_type=aggregation_type,
        )
        if property_scorers:
            stage.property_scorers.set(property_scorers)
        return stage

    def _ensure_net_checkpoint(self, net: models.ReinventNet):
        net.prepareData()
        ckpt_path = net.run_transfer_learning(device="cpu")
        self.assertTrue(os.path.isfile(ckpt_path))
        self.assertGreater(os.path.getsize(ckpt_path), 0)
        return ckpt_path

    # ----- TOML build tests -----

    def test_build_staged_toml_action_endpoint(self):
        net, _ = self._create_reinvent_net_via_api(build=False)
        self._ensure_net_checkpoint(net)

        env, agent, gen = self._mk_env_agent_gen(net, add_diversity=True)
        ps = self._create_property_scorer()
        self._create_stage_with_scorers(gen, max_steps=5, property_scorers=[ps])

        url = reverse("reinvent-build-toml", args=[gen.id])
        resp = self.client.post(url, data={"device": "cpu"}, format="json")
        self.assertEqual(resp.status_code, status.HTTP_201_CREATED, msg=resp.data)

        toml_path = resp.data["toml_path"]
        self.assertTrue(os.path.isfile(toml_path))

        cfg = open(toml_path, "r", encoding="utf-8").read()
        self.assertIn('run_type = "staged_learning"', cfg)
        self.assertIn("[parameters]", cfg)
        self.assertIn("[learning_strategy]", cfg)
        self.assertIn("[[stage]]", cfg)
        self.assertIn("[stage.scoring]", cfg)

    def test_build_staged_toml_includes_property_scorer(self):
        """Property scorer assigned to a stage appears in the TOML."""
        net, _ = self._create_reinvent_net_via_api(build=False)
        self._ensure_net_checkpoint(net)

        env, agent, gen = self._mk_env_agent_gen(net)
        ps = self._create_property_scorer(name="logP_test", property_name="SlogP", weight=0.8)
        self._create_stage_with_scorers(gen, property_scorers=[ps])

        toml_path = gen.build_staged_toml(device="cpu")
        cfg = open(toml_path, "r", encoding="utf-8").read()

        self.assertIn("SlogP", cfg)
        self.assertIn('name = "logP_test"', cfg)
        self.assertIn("weight = 0.8", cfg)

    def test_build_staged_toml_multiple_scorers_on_stage(self):
        """Multiple property scorers on one stage all appear in the TOML."""
        net, _ = self._create_reinvent_net_via_api(build=False)
        self._ensure_net_checkpoint(net)

        env, agent, gen = self._mk_env_agent_gen(net)
        ps1 = self._create_property_scorer(name="QED_scorer", property_name="Qed", weight=1.0)
        ps2 = self._create_property_scorer(name="MW_scorer", property_name="MolecularWeight", weight=0.5)
        self._create_stage_with_scorers(gen, property_scorers=[ps1, ps2])

        toml_path = gen.build_staged_toml(device="cpu")
        cfg = open(toml_path, "r", encoding="utf-8").read()

        self.assertIn("Qed", cfg)
        self.assertIn("MolecularWeight", cfg)
        self.assertIn('name = "QED_scorer"', cfg)
        self.assertIn('name = "MW_scorer"', cfg)

    def test_build_staged_toml_default_bad_smarts_only_in_first_stage(self):
        """Default bad SMARTS penalty is included only in stage 0."""
        net, _ = self._create_reinvent_net_via_api(build=False)
        self._ensure_net_checkpoint(net)

        env, agent, gen = self._mk_env_agent_gen(net)
        ps = self._create_property_scorer()
        self._create_stage_with_scorers(gen, order=0, property_scorers=[ps])
        self._create_stage_with_scorers(gen, order=1, max_steps=3, property_scorers=[ps])

        toml_path = gen.build_staged_toml(device="cpu")
        cfg = open(toml_path, "r", encoding="utf-8").read()

        # Should appear exactly once (first stage)
        count = cfg.count('name = "Unwanted SMARTS (default)"')
        self.assertEqual(count, 1, f"Expected 1 default bad SMARTS block, found {count}")

    def test_build_staged_toml_respects_bad_smarts_weight(self):
        """bad_smarts_weight from Reinvent model is reflected in TOML."""
        net, _ = self._create_reinvent_net_via_api(build=False)
        self._ensure_net_checkpoint(net)

        env, agent, gen = self._mk_env_agent_gen(net)
        gen.bad_smarts_weight = 0.7
        gen.save(update_fields=["bad_smarts_weight"])

        ps = self._create_property_scorer()
        self._create_stage_with_scorers(gen, property_scorers=[ps])

        toml_path = gen.build_staged_toml(device="cpu")
        cfg = open(toml_path, "r", encoding="utf-8").read()
        self.assertIn("weight = 0.7", cfg)

    def test_build_staged_toml_no_stages_raises(self):
        """build_staged_toml raises when no stages are defined."""
        net, _ = self._create_reinvent_net_via_api(build=False)
        self._ensure_net_checkpoint(net)

        env, agent, gen = self._mk_env_agent_gen(net)
        with self.assertRaises(RuntimeError):
            gen.build_staged_toml(device="cpu")

    def test_build_staged_toml_aggregation_type_per_stage(self):
        """Stage-level aggregation_type is used directly in the generated TOML."""
        net, _ = self._create_reinvent_net_via_api(build=False)
        self._ensure_net_checkpoint(net)

        env, agent, gen = self._mk_env_agent_gen(net)
        ps = self._create_property_scorer()
        self._create_stage_with_scorers(
            gen, property_scorers=[ps],
            aggregation_type="arithmetic_mean",
        )

        toml_path = gen.build_staged_toml(device="cpu")
        cfg = open(toml_path, "r", encoding="utf-8").read()
        self.assertIn('type = "arithmetic_mean"', cfg)

    # ----- RL run tests -----

    def test_run_staged_learning_via_celery_task_writes_rl_log(self):
        net, _ = self._create_reinvent_net_via_api(build=False)
        self._ensure_net_checkpoint(net)

        env, agent, gen = self._mk_env_agent_gen(net, add_diversity=False)
        ps = self._create_property_scorer()
        self._create_stage_with_scorers(gen, max_steps=3, property_scorers=[ps])

        # Execute the real task (eager mode => runs inline)
        res = runReinventStagedLearning.delay(gen.id, device="cpu")
        out = res.get()

        toml_path = out.get("toml_path")
        self.assertTrue(toml_path and os.path.isfile(toml_path))

        # RL log is written by ReinventAgent.run_staged_learning
        log_path = out.get("rl_log_path") or agent.get_rl_log_path()
        self.assertTrue(log_path and os.path.isfile(log_path))
        log_txt = open(log_path, "r", encoding="utf-8").read()
        self.assertIn("[CMD]", log_txt)
        self.assertTrue(len(log_txt.strip()) > 0)

    def test_run_staged_learning_endpoint_enqueues_task(self):
        net, _ = self._create_reinvent_net_via_api(build=False)
        self._ensure_net_checkpoint(net)

        env, agent, gen = self._mk_env_agent_gen(net, add_diversity=False)
        ps = self._create_property_scorer()
        self._create_stage_with_scorers(gen, max_steps=3, property_scorers=[ps])

        url = reverse("reinvent-run-staged-learning", args=[gen.id])
        resp = self.client.post(url, data={"device": "cpu"}, format="json")
        self.assertEqual(resp.status_code, status.HTTP_202_ACCEPTED, msg=resp.data)

        # In eager mode the task executes immediately; check RL log exists.
        log_path = agent.get_rl_log_path()
        self.assertTrue(os.path.isfile(log_path))
        txt = open(log_path, "r", encoding="utf-8").read()
        self.assertIn("[CMD]", txt)


# ---------------------------------------------------------------------
# Tests: get_csv_path picks latest stage CSV
# ---------------------------------------------------------------------
@override_settings(
    ROOT_URLCONF="genui.urls",
    CELERY_TASK_ALWAYS_EAGER=True,
    CELERY_TASK_EAGER_PROPAGATES=True,
)
class GetCsvPathTests(SetUpReinventMixIn, APITestCase):
    def test_get_csv_path_returns_highest_numbered_csv(self):
        """get_csv_path returns the CSV with the highest stage-number suffix."""
        net, _ = self._create_reinvent_net_via_api(build=False)
        self._ensure_net_checkpoint(net)

        env = self._create_dataset_like(
            models.ReinventEnvironment,
            name="CSV Path Test Env",
            prior_net=net,
            agent_net=net,
            aggregation_type="geometric_mean",
        )
        train_cfg = self._create_strategy_like(
            models.ReinventAgentTraining,
            model_instance=net,
            summary_csv_prefix="reinvent",
        )
        agent = self._create_model_like(
            models.ReinventAgent,
            name="CSV Path Agent",
            environment=env,
            training=train_cfg,
        )
        gen = self._create_model_like(
            models.Reinvent,
            name="CSV Path Run",
            environment=env,
            agent=agent,
        )

        sl_dir = agent._sl_dir()
        os.makedirs(sl_dir, exist_ok=True)

        prefix = agent._results_prefix(gen)
        # Write two fake CSVs: _1 (small) and _2 (larger)
        csv1 = os.path.join(sl_dir, f"{prefix}_1.csv")
        csv2 = os.path.join(sl_dir, f"{prefix}_2.csv")
        with open(csv1, "w") as f:
            f.write("step,SMILES,Score\n1,CCO,0.5\n")
        with open(csv2, "w") as f:
            f.write("step,SMILES,Score\n1,c1ccccc1,0.9\n2,CC(=O)O,0.8\n")

        result = agent.get_csv_path(generator=gen)
        self.assertEqual(result, csv2, f"Expected latest stage CSV ({csv2}), got {result}")

    def test_get_csv_path_fallback_plain_csv(self):
        """get_csv_path returns plain (unnumbered) CSV when no numbered files exist."""
        net, _ = self._create_reinvent_net_via_api(build=False)
        self._ensure_net_checkpoint(net)

        env = self._create_dataset_like(
            models.ReinventEnvironment,
            name="CSV Fallback Env",
            prior_net=net,
            agent_net=net,
            aggregation_type="geometric_mean",
        )
        train_cfg = self._create_strategy_like(
            models.ReinventAgentTraining,
            model_instance=net,
            summary_csv_prefix="reinvent",
        )
        agent = self._create_model_like(
            models.ReinventAgent,
            name="CSV Fallback Agent",
            environment=env,
            training=train_cfg,
        )
        gen = self._create_model_like(
            models.Reinvent,
            name="CSV Fallback Run",
            environment=env,
            agent=agent,
        )

        sl_dir = agent._sl_dir()
        os.makedirs(sl_dir, exist_ok=True)

        prefix = agent._results_prefix(gen)
        plain = os.path.join(sl_dir, f"{prefix}.csv")
        with open(plain, "w") as f:
            f.write("step,SMILES,Score\n1,CCO,0.5\n")

        result = agent.get_csv_path(generator=gen)
        self.assertEqual(result, plain)

    def _ensure_net_checkpoint(self, net):
        net.prepareData()
        net.run_transfer_learning(device="cpu")


# ---------------------------------------------------------------------
# Tests: Stage queryset filtering by generator
# ---------------------------------------------------------------------
@override_settings(
    ROOT_URLCONF="genui.urls",
    CELERY_TASK_ALWAYS_EAGER=True,
    CELERY_TASK_EAGER_PROPAGATES=True,
)
class StageFilteringTests(SetUpReinventMixIn, APITestCase):
    def test_stages_filtered_by_generator_param(self):
        """GET /reinvent/stages/?generator=X returns only stages for that run."""
        net, _ = self._create_reinvent_net_via_api(build=False)

        env = self._create_dataset_like(
            models.ReinventEnvironment,
            name="Stage Filter Env",
            prior_net=net,
            agent_net=net,
            aggregation_type="geometric_mean",
        )
        train_cfg = self._create_strategy_like(
            models.ReinventAgentTraining,
            model_instance=net,
        )
        agent = self._create_model_like(
            models.ReinventAgent,
            name="Stage Filter Agent",
            environment=env,
            training=train_cfg,
        )
        gen1 = self._create_model_like(
            models.Reinvent, name="Run A", environment=env, agent=agent,
        )
        gen2 = self._create_model_like(
            models.Reinvent, name="Run B", environment=env, agent=agent,
        )

        # Create 2 stages for gen1, 1 stage for gen2
        models.ReinventStage.objects.create(generator=gen1, order=0)
        models.ReinventStage.objects.create(generator=gen1, order=1)
        models.ReinventStage.objects.create(generator=gen2, order=0)

        url = reverse("reinventstage-list")
        resp1 = self.client.get(url, {"generator": gen1.id})
        self.assertEqual(resp1.status_code, status.HTTP_200_OK)
        self.assertEqual(len(resp1.data), 2)

        resp2 = self.client.get(url, {"generator": gen2.id})
        self.assertEqual(resp2.status_code, status.HTTP_200_OK)
        self.assertEqual(len(resp2.data), 1)


# ---------------------------------------------------------------------
# Tests: PropertyScorer API CRUD
# ---------------------------------------------------------------------
@override_settings(
    ROOT_URLCONF="genui.urls",
    CELERY_TASK_ALWAYS_EAGER=True,
    CELERY_TASK_EAGER_PROPAGATES=True,
)
class PropertyScorerAPITests(SetUpReinventMixIn, APITestCase):
    def test_create_property_scorer(self):
        url = reverse("reinvent-property-scorer-list")
        payload = {
            "name": "QED Scorer",
            "property_name": "Qed",
            "weight": 1.0,
            "project": self.project.id,
        }
        resp = self.client.post(url, data=payload, format="json")
        self.assertEqual(resp.status_code, status.HTTP_201_CREATED, msg=resp.data)
        self.assertEqual(resp.data["property_name"], "Qed")

    def test_list_property_scorers_filtered_by_project(self):
        models.PropertyScorer.objects.create(
            name="PS1", property_name="Qed", weight=1.0, project=self.project,
        )
        url = reverse("reinvent-property-scorer-list")
        resp = self.client.get(url, {"project_id": self.project.id})
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        self.assertGreaterEqual(len(resp.data), 1)

    def test_delete_property_scorer(self):
        ps = models.PropertyScorer.objects.create(
            name="ToDelete", property_name="SlogP", weight=0.5, project=self.project,
        )
        url = reverse("reinvent-property-scorer-detail", args=[ps.id])
        resp = self.client.delete(url)
        self.assertEqual(resp.status_code, status.HTTP_204_NO_CONTENT)
        self.assertFalse(models.PropertyScorer.objects.filter(pk=ps.id).exists())


