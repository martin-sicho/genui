# genui/generators/extensions/genuireinvent/models.py

from __future__ import annotations
import pkgutil
import inspect
import importlib

import os
import shutil
import subprocess
import tempfile
from typing import Tuple
import re
import csv

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import random

from django.conf import settings
from django.core.files.base import ContentFile
from django.db import models, close_old_connections
from django.utils import timezone

from genui.compounds.models import MolSet, ActivitySet
from genui.models.models import Model, ModelFile, TrainingStrategy, ValidationStrategy, ModelPerfomanceNN, ModelPerformance
from genui.projects.models import DataSet
from genui.generators.models import Generator


DEFAULT_PRIOR_REL = os.path.join("checkpoints", "prior", "reinvent.prior")

def _resolve_reinvent_prior_path() -> str:
    candidates = []
    candidates.append(getattr(settings, "REINVENT_PRIOR_PATH", None))
    try:
        candidates.append(getattr(settings, "GENUI_SETTINGS", {}).get("REINVENT_PRIOR_PATH"))
    except Exception:
        candidates.append(None)
    candidates.append(os.environ.get("REINVENT_PRIOR"))

    files_dir = None
    try:
        files_dir = getattr(settings, "GENUI_SETTINGS", {}).get("FILES_DIR")
    except Exception:
        files_dir = None
    if files_dir:
        candidates.append(os.path.join(files_dir, DEFAULT_PRIOR_REL))

    tried = [c for c in candidates if c]
    for c in tried:
        if os.path.isfile(c):
            return c

    raise FileNotFoundError(
        "REINVENT prior not found. Tried: "
        + ", ".join(tried or ["<no candidates>"])
        + ". Configure REINVENT_PRIOR (env) or REINVENT_PRIOR_PATH (settings), "
        + "or place the prior at "
        + (os.path.join(files_dir, DEFAULT_PRIOR_REL) if files_dir else "<FILES_DIR>/" + DEFAULT_PRIOR_REL)
    )

_BEST_EPOCH_RE = re.compile(
    r"Best\s+validation\s+loss\s*\(\s*(?P<loss>[-+]?(\d+(\.\d+)?|\.\d+))\s*\)\s*was\s*at\s*epoch\s*(?P<epoch>\d+)",
    re.IGNORECASE,
)

def _parse_best_from_log(text: str) -> tuple[int | None, float | None]:
    if not text:
        return None, None
    m = _BEST_EPOCH_RE.search(text)
    if not m:
        return None, None
    return int(m.group("epoch")), float(m.group("loss"))


# ───────────────────────────────────────────────────────────────────────────────
# Small helper: overwrite a hashed ModelFile in-place
# ───────────────────────────────────────────────────────────────────────────────
def _overwrite_filefield(mf: ModelFile, data: bytes | str, *, filename: str | None = None) -> None:
    """
    Overwrite an existing FileField content while keeping its hashed location.
    """
    if isinstance(data, str):
        data = data.encode("utf-8")

    # Prefer true in-place overwrite for local storage
    try:
        p = mf.file.path
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "wb") as fh:
            fh.write(data)
        return
    except Exception:
        pass

    # Fallback for storages without .path (S3 etc.)
    current_rel = mf.file.name
    if not current_rel:
        # first save
        mf.file.save(filename or f"aux_{mf.pk}", ContentFile(data), save=True)
        return

    try:
        mf.file.storage.delete(current_rel)
    except Exception:
        pass

    # Force same name
    mf.file.save(current_rel, ContentFile(data), save=False)
    mf.file.name = current_rel
    mf.save(update_fields=["file"])

def _bemis_murcko(smiles: str) -> str:
    m = Chem.MolFromSmiles(smiles)
    if not m: return ""
    core = MurckoScaffold.GetScaffoldForMol(m)
    return Chem.MolToSmiles(core, isomericSmiles=False) if core else ""

def _split_indices(n, frac, seed):
    r = random.Random(seed)
    idx = list(range(n))
    r.shuffle(idx)
    cut = max(1, int(n * frac))
    valid = set(idx[:cut])
    train = [i for i in idx if i not in valid]
    valid = list(valid)
    return train, valid


class _ReinventCLIModel:
    """
    Tiny façade so the builder/algorithm API works while training happens via CLI.
    """
    def __init__(self, net: "ReinventNet"):
        self._net = net
        self._checkpoint: str | None = None

    def fit(self, X=None, y=None):
        self._net.prepareData()
        self._checkpoint = self._net.run_transfer_learning(device="cpu")
        return self

    def loadStatesFromFile(self, path: str):
        return self

    def getModel(self):
        return {"checkpoint": self._checkpoint}


class ReinventNet(Model):
    # AUX notes (DrugEx-style)
    CORPUS_FULL_NOTE     = "reinvent_corpus_full"     # full cleaned .smi for CLI
    CORPUS_PREVIEW_NOTE  = "reinvent_corpus_preview"  # short preview for UI/tests
    TOML_FILE_NOTE       = "reinvent_tl_toml"         # generated TL config
    TRAIN_LOG_NOTE       = "reinvent_train_log"       # TL stdout/stderr log
    CHECKPOINT_FILE_NOTE = "reinvent_tl_checkpoint"   # where REINVENT writes
    CORPUS_TRAIN_NOTE = "reinvent_corpus_train"
    CORPUS_VALID_NOTE = "reinvent_corpus_valid"
    PRIOR_FILE_NOTE = "reinvent_prior_copy"

    molset = models.ForeignKey(MolSet, on_delete=models.CASCADE, null=True)
    parent = models.ForeignKey("self", on_delete=models.CASCADE, null=True)

    # ── AUX getters (create the record lazily with empty payload) ──────────────
    def _get_or_create_aux(self, note: str, filename: str) -> ModelFile:
        mf = self.files.filter(kind=ModelFile.AUXILIARY, note=note).first()
        if mf is None:
            mf = ModelFile.create(self, filename, ContentFile(b""), note=note)
        return mf

    @property
    def corpusFileTrain(self):  # backwards-compat alias
        return self.corpusTrainFile

    @property
    def corpusTrainFile(self) -> ModelFile:
        return self._get_or_create_aux(self.CORPUS_TRAIN_NOTE, f"corpus_train_{self.pk}.smi")

    @property
    def corpusValidFile(self) -> ModelFile:
        return self._get_or_create_aux(self.CORPUS_VALID_NOTE, f"corpus_valid_{self.pk}.smi")

    @property
    def corpusFullFile(self) -> ModelFile:
        # Full cleaned corpus consumed by REINVENT CLI
        return self._get_or_create_aux(self.CORPUS_FULL_NOTE, f"corpus_full_{self.pk}.smi")

    @property
    def corpusPreviewFile(self) -> ModelFile:
        # Optional short preview for UI/tests
        return self._get_or_create_aux(self.CORPUS_PREVIEW_NOTE, f"corpus_preview_{self.pk}.smi")

    @property
    def tlTomlFile(self) -> ModelFile:
        return self._get_or_create_aux(self.TOML_FILE_NOTE, f"tl_reinvent_{self.pk}.toml")

    @property
    def trainLogFile(self) -> ModelFile:
        return self._get_or_create_aux(self.TRAIN_LOG_NOTE, f"reinvent_training_{self.pk}.log")

    @property
    def checkpointFile(self) -> ModelFile:
        # We keep the checkpoint managed as an AUX file too
        return self._get_or_create_aux(self.CHECKPOINT_FILE_NOTE, f"reinvent_{self.pk}.model")

    # Backwards-compat convenience (tests may call this):
    def get_clean_corpus_path(self) -> str:
        return self.corpusFullFile.path

    # ── Prior path  ────────────────────────────────────────────────
    @property
    def priorFile(self) -> ModelFile:
        # stored as an AUX file tied to this model
        return self._get_or_create_aux(self.PRIOR_FILE_NOTE, f"reinvent_prior_{self.pk}.prior")

    def ensure_prior_copy(self) -> ModelFile:
        """
        Make sure this model has its own prior copy (for reproducibility).
        """
        src = _resolve_reinvent_prior_path()
        with open(src, "rb") as f:
            data = f.read()
        _overwrite_filefield(self.priorFile, data, filename=os.path.basename(src))
        return self.priorFile

    def get_prior_path(self) -> str:
        """
        Prefer the project/model-owned copy if present, else fall back to global.
        """
        mf = self.files.filter(kind=ModelFile.AUXILIARY, note=self.PRIOR_FILE_NOTE).first()
        if mf and mf.file:
            try:
                # local storage
                return mf.file.path
            except Exception:
                pass
        # fallback (global location ensured by genuisetup)
        return _resolve_reinvent_prior_path()

    # ── Clean corpus preparation (hashed AUX only) ─────────────────────────────
    def prepareData(self) -> Tuple[ModelFile, ModelFile]:
        """
        Clean SMILES via reinvent.datapipeline and write:
          - Full cleaned corpus directly to corpusFullFile.path (hashed in media/)
          - Short preview (first 1000 lines) into corpusPreviewFile (hashed)
        """
        if not self.molset:
            raise RuntimeError(f"No MolSet attached to {self}.")

        # Decide input for datapipeline
        input_path = None
        if getattr(self.molset, "files", None) and self.molset.files.exists():
            f = self.molset.files.first()
            if f and getattr(f, "file", None):
                input_path = f.file.path

        # If needed, emit a temporary TSV with a SMILES header
        temp_in = None
        if not input_path:
            with tempfile.NamedTemporaryFile(prefix=f"reinvent_raw_{self.pk}_", suffix=".smi.tsv", delete=False) as tf:
                temp_in = tf.name
            with open(temp_in, "w", encoding="utf-8") as w:
                w.write("SMILES\n")
                for s in self.molset.allSmiles:
                    w.write(s + "\n")
            input_path = temp_in

        out_full_path = self.corpusFullFile.path  # hashed media path

        try:
            from reinvent.datapipeline import preprocess
        except Exception as e:
            if temp_in:
                try:
                    os.remove(temp_in)
                except OSError:
                    pass
            raise RuntimeError("reinvent.datapipeline.preprocess is required.") from e

        cfg_text = f"""\
        input_csv_file = "{input_path}"
        smiles_column = "SMILES"
        separator = "\\t"
        output_smiles_file = "{out_full_path}"

        [filter]
        elements = []
        transforms = ["standard"]
        inchi_key_deduplicate = true
        """
        with tempfile.NamedTemporaryFile(prefix=f"reinvent_preprocess_{self.pk}_",
                                         suffix=".toml", delete=False) as tf:
            cfg_path = tf.name
        try:
            with open(cfg_path, "w", encoding="utf-8") as fh:
                fh.write(cfg_text)
            args = type("Args", (), {"config_filename": cfg_path, "log_filename": None})
            preprocess.main(args)
        finally:
            try:
                os.remove(cfg_path)
            except OSError:
                pass
            if temp_in:
                try:
                    os.remove(temp_in)
                except OSError:
                    pass

        # 2) Read CLEANED full corpus and split
        vs = getattr(self, "validationStrategy", None)
        method = (getattr(vs, "split_method", None) or "random").lower()
        frac = max(0.0, min(0.9, float(getattr(vs, "valid_fraction", 0.1))))
        seed = int(getattr(vs, "random_seed", 1337))
        cutoff = getattr(vs, "temporal_cutoff", None)
        max_valid = int(getattr(vs, "validSetSize", 0)) or None

        with open(out_full_path, "r", encoding="utf-8") as fh:
            smiles = [ln.strip() for ln in fh if ln.strip()]
        # Guard empty corpus before splitting. If the preprocessor yields 0–1 lines, your split can produce empty files.
        if not smiles:
            raise RuntimeError(f"Cleaned corpus is empty at {out_full_path}.")
        if len(smiles) == 1:
            _overwrite_filefield(self.corpusTrainFile, smiles[0] + "\n",
                                 filename=os.path.basename(self.corpusTrainFile.file.name))
            _overwrite_filefield(self.corpusValidFile, "",
                                 filename=os.path.basename(self.corpusValidFile.file.name))
            # preview build as you do…
            return self.corpusTrainFile, self.corpusValidFile

        if method == "scaffold":
            buckets = {}
            for s in smiles:
                scf = _bemis_murcko(s) or f"NOSCAF_{hash(s) % 10_000_000}"
                buckets.setdefault(scf, []).append(s)
            rng = random.Random(seed)
            scaf_ids = list(buckets.keys());
            rng.shuffle(scaf_ids)
            valid_target = max(1, int(len(smiles) * frac))
            train, valid, acc = [], [], 0
            for scf in scaf_ids:
                grp = buckets[scf]
                if acc < valid_target:
                    valid.extend(grp);
                    acc += len(grp)
                else:
                    train.extend(grp)
        elif method == "temporal" and cutoff:
            raise NotImplementedError("Temporal split needs SMILES->date mapping in MolSet.")
        else:
            tr_idx, va_idx = _split_indices(len(smiles), frac, seed)
            train = [smiles[i] for i in tr_idx]
            valid = [smiles[i] for i in va_idx]

        if max_valid is not None and len(valid) > max_valid:
            valid = valid[:max_valid]
        if not train:
            move_n = max(1, len(valid) // 2)
            train, valid = valid[:move_n], valid[move_n:]

        _overwrite_filefield(self.corpusTrainFile, "\n".join(train) + "\n",
                             filename=os.path.basename(self.corpusTrainFile.file.name))
        _overwrite_filefield(self.corpusValidFile, "\n".join(valid) + "\n",
                             filename=os.path.basename(self.corpusValidFile.file.name))

        # 3) Build preview from CLEANED corpus
        head = []
        with open(out_full_path, "r", encoding="utf-8") as f:
            for i, ln in enumerate(f):
                if i >= 1000: break
                s = ln.strip()
                if s: head.append(s)
        preview_text = ("\n".join(head) + "\n") if head else ""
        _overwrite_filefield(self.corpusPreviewFile, preview_text,
                             filename=os.path.basename(self.corpusPreviewFile.file.name))

        # 4) Return actual train/valid
        return self.corpusTrainFile, self.corpusValidFile

    # ── TOML (hashed AUX only) ────────────────────────────────────────────────
    def build_tl_toml(self, *, device: str = "cpu") -> str:
        ts = self.trainingStrategy
        if not isinstance(ts, ReinventNetTraining):
            raise RuntimeError("ReinventNetTraining required.")

        prior = self.get_prior_path()
        out_path = self.checkpointFile.path   # REINVENT will write here
        train = self.corpusTrainFile.path
        valid = self.corpusValidFile.path
        sbs = max(100, ts.sample_batch_size)
        tb_dir = os.path.join(settings.MEDIA_ROOT, "models", f"tb_TL_{self.pk}")
        os.makedirs(tb_dir, exist_ok=True)

        body = f"""\
run_type = "transfer_learning"
device = "{device}"
tb_logdir = "{tb_dir}"

[parameters]
num_epochs = {ts.epochs}
save_every_n_epochs = {ts.save_every_n_epochs}
batch_size = {ts.batch_size}
sample_batch_size = {sbs}

input_model_file = "{prior}"
smiles_file = "{train}"
validation_smiles_file = "{valid}"
output_model_file = "{out_path}"
"""

        _overwrite_filefield(
            self.tlTomlFile,
            body,
            filename=os.path.basename(self.tlTomlFile.file.name),
        )
        return self.tlTomlFile.path

    @staticmethod
    def _pick_best_checkpoint(tb_dir: str) -> tuple[int, float] | None:
        try:
            from tensorboard.backend.event_processing.event_accumulator import EventAccumulator
            ea = EventAccumulator(tb_dir);
            ea.Reload()
            vals = ea.Scalars("valid/nll") or ea.Scalars("validation/nll")
            if not vals: return None
            best = min(vals, key=lambda x: x.value)
            return (best.step, best.value)
        except Exception:
            return None

    def get_active_checkpoint_path(self) -> str:
        """
        Returns the canonical checkpoint to load for the next stage.
        Prefer the selected best-epoch (copied into checkpointFile.path).
        """
        p = self.checkpointFile.path
        if not os.path.isfile(p):
            raise FileNotFoundError(f"Active checkpoint missing at {p}.")
        return p

    # ── TL run (hashed AUX only) ───────────────────────────────────────────────
    def run_transfer_learning(self, *, device: str = "cpu") -> str:
        toml_path = self.build_tl_toml(device=device)
        out_path = self.checkpointFile.path

        reinvent_bin = (getattr(settings, "REINVENT_BIN", None)
                        or os.environ.get("REINVENT_BIN")
                        or shutil.which("reinvent"))
        if not reinvent_bin:
            raise RuntimeError("REINVENT binary not found. Set settings.REINVENT_BIN or $REINVENT_BIN.")

        cmd = [reinvent_bin, toml_path]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, bufsize=1, cwd=settings.BASE_DIR)
        lines = [ln for ln in (proc.stdout or [])]
        rc = proc.wait()

        log_text = f"[CMD] {' '.join(cmd)}\n{''.join(lines)}"
        _overwrite_filefield(self.trainLogFile, log_text,
                             filename=os.path.basename(self.trainLogFile.file.name))

        if rc != 0:
            raise RuntimeError(f"REINVENT TL failed (exit={rc}). See TOML: {toml_path}")

        # optional: swap to best epoch checkpoint based on TensorBoard
        tb_dir = os.path.join(settings.MEDIA_ROOT, "models", f"tb_TL_{self.pk}")
        best_epoch, best_loss = _parse_best_from_log(log_text)

        if best_epoch is not None and best_loss is not None:
            ts = self.trainingStrategy
            ts.best_epoch = best_epoch
            ts.best_valid_loss = best_loss
            ts.save(update_fields=["best_epoch", "best_valid_loss"])

        return out_path

    # Keep the façade so builders can call into “a model”
    def getModel(self):
        return _ReinventCLIModel(self)


class ReinventNetValidation(ValidationStrategy):
    validSetSize = models.IntegerField(default=10000)  # keep if you want “cap”
    split_method = models.CharField(
        max_length=16, default="random",  # "random" | "scaffold" | "temporal"
    )
    valid_fraction = models.FloatField(default=0.1)  # ignored if validSetSize used
    random_seed = models.IntegerField(default=1337)
    temporal_cutoff = models.CharField(max_length=32, null=True, blank=True)  # e.g. "2024-06-01"


class ReinventNetTraining(TrainingStrategy):
    epochs = models.IntegerField(default=10)
    batch_size = models.IntegerField(default=64)
    save_every_n_epochs = models.IntegerField(default=1)
    sample_batch_size = models.IntegerField(default=100)

    best_epoch = models.IntegerField(null=True, blank=True)
    best_valid_loss = models.FloatField(null=True, blank=True)

    def processMetaData(self, metadata: dict):
            self.epochs = metadata.get("epochs", self.epochs)
            self.batch_size = metadata.get("batch_size", self.batch_size)
            self.sample_batch_size = metadata.get("sample_batch_size", self.sample_batch_size)
            self.save()


# =====================================================================
#  STAGED LEARNING MODULE  (AFTER TL CODE)
# =====================================================================

# Hardcoded unwanted SMARTS from REINVENT supplement
UNWANTED_SMARTS_DEFAULT = [
    "[*;r8]", "[*;r9]", "[*;r10]", "[*;r11]", "[*;r12]", "[*;r13]", "[*;r14]",
    "[*;r15]", "[*;r16]", "[*;r17]",
    "[#8][#8]", "[#6;+]", "[#16][#16]",
    "[#7;!n][S;!$(S(=O)=O)]", "[#7;!n][#7;!n]",
    "C#C", "C(=[O,S])[O,S]",
    "[#7;!n][C;!$(C(=[O,N])[N,O])][#16;!s]",
    "[#7;!n][C;!$(C(=[O,N])[N,O])][#7;!n]",
    "[#7;!n][C;!$(C(=[O,N])[N,O])][#8;!o]",
    "[#8;!o][C;!$(C(=[O,N])[N,O])][#16;!s]",
    "[#8;!o][C;!$(C(=[O,N])[N,O])][#8;!o]",
    "[#16;!s][C;!$(C(=[O,N])[N,O])][#16;!s]",
]


class ReinventEnvironmentHelper:

    @staticmethod
    def get_diversity_filters():
        try:
            import reinvent.runmodes.RL.memories as mem
            from reinvent.runmodes.RL.memories.diversity_filter import DiversityFilter
        except ModuleNotFoundError:
            return []
        results = []

        for _, modname, _ in pkgutil.walk_packages(mem.__path__, prefix="reinvent.runmodes.RL.memories."):
            lname = modname.lower()
            if any(x in lname for x in ["murcko", "topological", "similarity", "penalize"]):
                module = importlib.import_module(modname)
                for name, obj in inspect.getmembers(module, inspect.isclass):
                    if issubclass(obj, DiversityFilter) and obj is not DiversityFilter:
                        results.append(name)

        return sorted(set(results))

    @staticmethod
    def get_learning_strategies():
        return ["dap"]


def df_choices():
    return [(x, x) for x in ReinventEnvironmentHelper.get_diversity_filters()]


# =====================================================================
# ENVIRONMENT
# =====================================================================

class ReinventDiversityFilter(models.Model):
    type = models.CharField(max_length=128, choices=df_choices)
    bucket_size = models.IntegerField(default=25)
    minscore = models.FloatField(default=0.4)
    minsimilarity = models.FloatField(default=0.4)
    penalty_multiplier = models.FloatField(default=0.5)

    def to_reinvent(self):
        d = {"type": self.type, "bucket_size": self.bucket_size, "minscore": self.minscore}
        if self.type == "ScaffoldSimilarity":
            d["minsimilarity"] = self.minsimilarity
        if self.type == "PenalizeSameSmiles":
            d["penalty_multiplier"] = self.penalty_multiplier
        return d

class ReinventEnvironment(DataSet):
    name = models.CharField(max_length=255)

    # Backwards-compatible (already in your model)
    prior_model = models.ForeignKey(
        ModelFile, on_delete=models.PROTECT, related_name="reinvent_prior_files",
        null=True, blank=True
    )
    agent_model = models.ForeignKey(
        ModelFile, on_delete=models.PROTECT, related_name="reinvent_agent_files",
        null=True, blank=True
    )

    # NEW: choose from existing models (user-added or previously trained)
    prior_net = models.ForeignKey(
        "ReinventNet", on_delete=models.PROTECT, related_name="reinvent_as_prior",
        null=True, blank=True
    )
    agent_net = models.ForeignKey(
        "ReinventNet", on_delete=models.PROTECT, related_name="reinvent_as_agent",
        null=True, blank=True
    )

    # Reward scheme aggregation type (directly stored, no separate model needed)
    aggregation_type = models.CharField(
        max_length=64,
        choices=[
            ("geometric_mean", "Geometric Mean (Balanced - all scores must be good)"),
            ("arithmetic_mean", "Arithmetic Mean (Flexible - based on weights)"),
        ],
        default="geometric_mean"
    )

    diversity_filter = models.ForeignKey(ReinventDiversityFilter, null=True, blank=True, on_delete=models.SET_NULL)

    inception_smiles = models.ForeignKey(ModelFile, null=True, blank=True, on_delete=models.SET_NULL)
    inception_memory_size = models.IntegerField(default=0)
    inception_sample_size = models.IntegerField(default=0)

    def get_prior_path(self) -> str:
        if self.prior_net:
            return self.prior_net.get_active_checkpoint_path()
        if self.prior_model:
            return self.prior_model.file.path
        raise RuntimeError("No prior selected: set prior_net or prior_model.")

    def get_agent_path(self) -> str:
        if self.agent_net:
            return self.agent_net.get_active_checkpoint_path()
        if self.agent_model:
            return self.agent_model.file.path
        raise RuntimeError("No agent selected: set agent_net or agent_model.")


# =====================================================================
# REINVENT4 TRANSFORM UTILITIES
# =====================================================================

# Default transforms keyed by REINVENT4 property component name.
# These are auto-applied when no custom transform is specified.
REINVENT4_TRANSFORM_DEFAULTS = {
    # ── QED (0-1, want high) ─────────────────────────────────────────────
    "Qed":               {"type": "sigmoid",        "low": 0.5,  "high": 0.9,  "k": 0.5},
    # ── Lipophilicity (want 0-5, optimal ~1-4) ───────────────────────────
    # window=5 → coef_si = 200/(0.15*5) ≈ 25 for smooth transition
    "SlogP":             {"type": "double_sigmoid",  "low": 0.0,  "high": 5.0,  "coef_div": 100.0, "coef_si": 25.0, "coef_se": 25.0},
    # ── Molecular weight (want 150-500) ──────────────────────────────────
    # window=350 → coef_si = 200/(0.15*350) ≈ 4 for smooth transition
    "MolecularWeight":   {"type": "double_sigmoid",  "low": 150.0,"high": 500.0,"coef_div": 100.0, "coef_si": 4.0,  "coef_se": 4.0},
    # ── TPSA (want < 140) ────────────────────────────────────────────────
    "TPSA":              {"type": "reverse_sigmoid", "low": 0.0,  "high": 140.0,"k": 0.5},
    # ── H-bond acceptors (want ≤ 10) ─────────────────────────────────────
    "HBondAcceptors":    {"type": "reverse_sigmoid", "low": 0.0,  "high": 10.0, "k": 0.5},
    # ── H-bond donors (want ≤ 5) ─────────────────────────────────────────
    "HBondDonors":       {"type": "reverse_sigmoid", "low": 0.0,  "high": 5.0,  "k": 0.5},
    # ── Rotatable bonds (want ≤ 10) ──────────────────────────────────────
    "NumRotBond":        {"type": "reverse_sigmoid", "low": 0.0,  "high": 10.0, "k": 0.5},
    # ── Csp3 (want high, > 0.3) ──────────────────────────────────────────
    "Csp3":              {"type": "sigmoid",         "low": 0.2,  "high": 0.8,  "k": 0.5},
    # ── Ring counts (want 1-4) ────────────────────────────────────────────
    # window=3 → coef_si = 200/(0.15*3) ≈ 45 for smooth transition
    "NumRings":          {"type": "double_sigmoid",  "low": 1.0,  "high": 4.0,  "coef_div": 100.0, "coef_si": 45.0, "coef_se": 45.0},
    "NumAromaticRings":  {"type": "double_sigmoid",  "low": 0.0,  "high": 3.0,  "coef_div": 100.0, "coef_si": 45.0, "coef_se": 45.0},
    "LargestRingSize":   {"type": "reverse_sigmoid", "low": 3.0,  "high": 8.0,  "k": 0.5},
    # ── SAScore (1–10, lower is better) ──────────────────────────────────
    "SAScore":           {"type": "reverse_sigmoid", "low": 1.0,  "high": 6.0,  "k": 0.5},
    # ── Tanimoto distance/similarity (want > 0.4) ────────────────────────
    "TanimotoSimilarity":{"type": "sigmoid",         "low": 0.3,  "high": 0.8,  "k": 0.5},
    "TanimotoDistance":  {"type": "sigmoid",         "low": 0.3,  "high": 0.8,  "k": 0.5},
    # ── Substructure (0 or 1 binary) ─────────────────────────────────────
    "MatchingSubstructure": {"type": "right_step",   "low": 0.5,  "high": 0.5},
    "GroupCount":        {"type": "double_sigmoid",  "low": 1.0,  "high": 5.0,  "coef_div": 100.0, "coef_si": 45.0, "coef_se": 45.0},
    # ── Atom counts ──────────────────────────────────────────────────────
    # window=20 → coef_si = 200/(0.15*20) ≈ 67 for ~15% transition
    "NumHeavyAtoms":     {"type": "double_sigmoid",  "low": 10.0, "high": 40.0, "coef_div": 100.0, "coef_si": 6.7, "coef_se": 6.7},
    # window=5 → coef_si ≈ 25
    "NumHeteroAtoms":    {"type": "double_sigmoid",  "low": 1.0,  "high": 6.0,  "coef_div": 100.0, "coef_si": 45.0, "coef_se": 45.0},
    # ── Aliphatic rings (want 0-3) ───────────────────────────────────────
    "NumAliphaticRings": {"type": "double_sigmoid",  "low": 0.0,  "high": 3.0,  "coef_div": 100.0, "coef_si": 45.0, "coef_se": 45.0},
    # ── Stereocenters (want 0-3, fewer is easier to synthesise) ─────────
    "NumAtomSteroCenters": {"type": "reverse_sigmoid", "low": 0.0, "high": 3.0, "k": 0.5},
}

# Human-readable metadata for each transform type (for API / UI)
REINVENT4_TRANSFORM_META = {
    "sigmoid": {
        "label": "Sigmoid",
        "description": "Smooth S-curve: reward rises from 0→1 as value passes from `low` to `high`.",
        "params": {
            "low":  {"type": "number", "title": "Low (x where output ≈ 0.05)", "default": 0.0},
            "high": {"type": "number", "title": "High (x where output ≈ 0.95)", "default": 1.0},
            "k":    {"type": "number", "title": "Steepness k", "default": 0.5, "min": 0.01, "max": 10.0},
        },
    },
    "reverse_sigmoid": {
        "label": "Reverse Sigmoid",
        "description": "Smooth S-curve inverted: reward falls from 1→0 as value passes from `low` to `high`.",
        "params": {
            "low":  {"type": "number", "title": "Low (x where output ≈ 0.95)", "default": 0.0},
            "high": {"type": "number", "title": "High (x where output ≈ 0.05)", "default": 1.0},
            "k":    {"type": "number", "title": "Steepness k", "default": 0.5, "min": 0.01, "max": 10.0},
        },
    },
    "double_sigmoid": {
        "label": "Double Sigmoid (Hump)",
        "description": "Bell-shaped curve: reward peaks between `low` and `high`, falls off on both sides.",
        "params": {
            "low":      {"type": "number", "title": "Low boundary", "default": 0.0},
            "high":     {"type": "number", "title": "High boundary", "default": 1.0},
            "coef_div": {"type": "number", "title": "coef_div (scale factor)", "default": 100.0},
            "coef_si":  {"type": "number", "title": "coef_si (left steepness)", "default": 150.0},
            "coef_se":  {"type": "number", "title": "coef_se (right steepness)", "default": 150.0},
        },
    },
    "right_step": {
        "label": "Right Step",
        "description": "Returns 1.0 for values ≥ `high`, 0.0 otherwise.",
        "params": {
            "low":  {"type": "number", "title": "Low (unused)", "default": 0.0},
            "high": {"type": "number", "title": "Threshold (step point)", "default": 0.5},
        },
    },
    "left_step": {
        "label": "Left Step",
        "description": "Returns 1.0 for values ≤ `low`, 0.0 otherwise.",
        "params": {
            "low":  {"type": "number", "title": "Threshold (step point)", "default": 0.5},
            "high": {"type": "number", "title": "High (unused)", "default": 1.0},
        },
    },
    "step": {
        "label": "Step (Window)",
        "description": "Returns 1.0 for values in [low, high], 0.0 outside.",
        "params": {
            "low":  {"type": "number", "title": "Low bound", "default": 0.0},
            "high": {"type": "number", "title": "High bound", "default": 1.0},
        },
    },
    "exponential_decay": {
        "label": "Exponential Decay",
        "description": "exp(-k·x) for x≥0, clipped to 1.0 for x<0. Penalises large positive values.",
        "params": {
            "k": {"type": "number", "title": "Decay rate k", "default": 1.0, "min": 0.001},
        },
    },
}


def compute_transform_values(transform_dict: dict, x_values: list) -> list:
    """
    Apply a REINVENT4 transform to a list of x values.
    Returns list of float outputs in [0, 1].
    """
    import numpy as np

    t = dict(transform_dict)
    t_type = t.get("type", "sigmoid").lower().replace("_", "")

    x = np.array(x_values, dtype=np.float64)

    if t_type == "sigmoid":
        low = float(t.get("low", 0.0))
        high = float(t.get("high", 1.0))
        k = float(t.get("k", 0.5))
        center = (high + low) / 2.0
        xc = x - center
        if (high - low) == 0:
            k_eff = 10.0 * k
            y = (k_eff * xc > 0).astype(np.float64)
        else:
            k_eff = 10.0 * k / (high - low)
            h = k_eff * xc * np.log(10)
            y = np.where(h >= 0,
                         1.0 / (1.0 + np.exp(-h)),
                         np.exp(h) / (1.0 + np.exp(h)))
        return y.tolist()

    elif t_type == "reversesigmoid":
        low = float(t.get("low", 0.0))
        high = float(t.get("high", 1.0))
        k = float(t.get("k", 0.5))
        center = (high + low) / 2.0
        xc = x - center
        if (high - low) == 0:
            k_eff = 10.0 * k
            y = (k_eff * xc > 0).astype(np.float64)
        else:
            k_eff = 10.0 * k / (high - low)
            h = k_eff * xc * np.log(10)
            y = np.where(h >= 0,
                         1.0 / (1.0 + np.exp(-h)),
                         np.exp(h) / (1.0 + np.exp(h)))
        return (1.0 - y).tolist()

    elif t_type == "doublesigmoid":
        from reinvent.scoring.transforms.sigmoid_functions import double_sigmoid
        low = float(t.get("low", 0.0))
        high = float(t.get("high", 1.0))
        coef_div = float(t.get("coef_div", 100.0))
        # Default coef_si/coef_se to a window-proportional value so preview is never invisible.
        # Transition width ≈ 2*coef_div/coef_si; we want ~15% of window by default.
        window = abs(high - low) or 1.0
        auto_coef = max(0.5, round(2.0 * coef_div / (0.15 * window), 1))
        coef_si = float(t.get("coef_si", auto_coef))
        coef_se = float(t.get("coef_se", auto_coef))
        # double_sigmoid(x, x_left, x_right, k, k_left, k_right)
        # where x_left=low, x_right=high, k=coef_div, k_left=coef_si, k_right=coef_se
        y = double_sigmoid(x, low, high, coef_div, coef_si, coef_se)
        return y.tolist()

    elif t_type == "rightstep":
        high = float(t.get("high", 0.5))
        return [1.0 if v >= high else 0.0 for v in x_values]

    elif t_type == "leftstep":
        low = float(t.get("low", 0.5))
        return [1.0 if v <= low else 0.0 for v in x_values]

    elif t_type == "step":
        low = float(t.get("low", 0.0))
        high = float(t.get("high", 1.0))
        return [1.0 if low <= v <= high else 0.0 for v in x_values]

    elif t_type == "exponentialdecay":
        k = float(t.get("k", 1.0))
        y = np.where(x < 0, 1.0, np.exp(-k * x))
        return y.tolist()

    # Fallback: identity
    return x_values


# =====================================================================
# SCORING
# =====================================================================

class ScoreModifier(DataSet):
    """
    Base class for score modifiers (DrugEx style).
    Stored polymorphically; concrete implementations are subclasses.
    """

    def to_reinvent_transform(self) -> dict:
        raise NotImplementedError("Override in subclass.")



class ClippedScore(ScoreModifier):
    """
    Mirrors DrugEx.modifiers.ClippedScore / SmoothClippedScore configuration,
    but exports REINVENT transform dict.
    """
    upper = models.FloatField(null=False)
    lower = models.FloatField(null=False, default=0.0)
    high = models.FloatField(null=False, default=1.0)
    low = models.FloatField(null=False, default=0.0)
    smooth = models.BooleanField(null=False, default=False)

    def to_reinvent_transform(self) -> dict:
        if not self.smooth:
            # Hard clipped effect via double sigmoid window
            return {
                "type": "double_sigmoid",
                "low": float(self.lower),
                "high": float(self.upper),
                "coef_div": float(self.upper - self.lower) if self.upper != self.lower else 1.0,
                "coef_si": 20,
                "coef_se": 20,
            }

        # smooth=True → REINVENT reverse_sigmoid
        # Semantics: raw scores near `high` → reward ≈ 1,
        #            raw scores near `low`  → reward ≈ 0.
        # k controls steepness of the sigmoid curve (0.1 = gentle, 1.0 = steep).
        # We derive k from the `high` and `low` output fields:
        #   high/low ∈ [0, 1] in the UI.  Use (1 - high + low) clamped to
        #   a sensible range as the "softness" knob: when high=1,low=0
        #   → k=0.25 (smooth); when high=0.9,low=0.1 → k=0.45 (steeper).
        k = max(0.05, (1.0 - float(self.high) + float(self.low)) + 0.25)
        return {
            "type": "reverse_sigmoid",
            "low": float(self.lower),
            "high": float(self.upper),
            "k": round(k, 6),
        }


class SmoothHump(ScoreModifier):
    """
    Mirrors DrugEx.modifiers.SmoothHump, exports REINVENT hump transform.
    """
    upper = models.FloatField(null=False, default=1.0)
    lower = models.FloatField(null=False, default=0.0)
    sigma = models.FloatField(null=False, default=0.5)

    def to_reinvent_transform(self) -> dict:
        coef_div = float(self.upper - self.lower) if self.upper != self.lower else 1.0
        return {
            "type": "double_sigmoid",
            "low": float(self.lower),
            "high": float(self.upper),
            "coef_div": coef_div,
            "coef_si": int(float(self.sigma) * 20),
            "coef_se": int(float(self.sigma) * 20),
        }


class ScoringMethod(models.Model):
    project = models.ForeignKey(
        "projects.Project", on_delete=models.CASCADE,
        null=True, blank=True, related_name="+"
    )
    name = models.CharField(max_length=255)
    weight = models.FloatField(default=1.0)
    modifier = models.ForeignKey("ScoreModifier", null=True, blank=True, on_delete=models.SET_NULL)

    class Meta:
        abstract = True

    def build_transform(self) -> dict | None:
        mod = getattr(self, "modifier", None)
        if not mod:
            return None
        # If ScoreModifier is base, try to downcast via reverse relations.
        # Polymorphic reverse OneToOne lookups raise DoesNotExist (not
        # AttributeError) when the child row doesn't exist, so we must
        # catch that explicitly.
        for attr in ("clippedscore", "smoothhump"):
            try:
                child = getattr(mod, attr, None)
            except Exception:
                child = None
            if child is not None:
                mod = child
                break
        fn = getattr(mod, "to_reinvent_transform", None)
        return fn() if callable(fn) else None


class PropertyScorer(ScoringMethod):
    property_name = models.CharField(max_length=64)
    # Component-level params for REINVENT 4 components that require them.
    # Stored as JSON, e.g. {"smiles": ["CCO"], "radius": 3, "use_counts": true}
    component_params = models.JSONField(null=True, blank=True, default=None)
    # REINVENT4 transform configuration stored as JSON.
    # e.g. {"type": "sigmoid", "low": 0.0, "high": 1.0, "k": 0.5}
    # When null, a sensible default transform is auto-applied per property.
    transform_params = models.JSONField(null=True, blank=True, default=None)

    def build_transform(self) -> dict | None:
        """
        Return transform dict for this scorer.
        Priority:
          1. Explicit transform_params stored on this instance
          2. Legacy modifier FK (ClippedScore / SmoothHump)
          3. Auto-default from REINVENT4_TRANSFORM_DEFAULTS keyed by property_name
        """
        # 1) Explicit transform_params
        if self.transform_params and isinstance(self.transform_params, dict):
            return dict(self.transform_params)

        # 2) Legacy modifier
        mod = getattr(self, "modifier", None)
        if mod:
            for attr in ("clippedscore", "smoothhump"):
                try:
                    child = getattr(mod, attr, None)
                except Exception:
                    child = None
                if child is not None:
                    mod = child
                    break
            fn = getattr(mod, "to_reinvent_transform", None)
            if callable(fn):
                return fn()

        # 3) Auto-default by property name
        return dict(REINVENT4_TRANSFORM_DEFAULTS.get(self.property_name, {})) or None

class GenUIModelScorer(ScoringMethod):
    model = models.ForeignKey(Model, on_delete=models.PROTECT)

class UnwantedSmartsScorer(ScoringMethod):
    enabled = models.BooleanField(default=True)

    def load_patterns(self):
        """
        Return the SMARTS patterns to use for custom alerts.

        For now we just use a static default list (UNWANTED_SMARTS_DEFAULT).
        This can be extended later to load user-defined SMARTS from another
        model or file if needed.
        """
        return UNWANTED_SMARTS_DEFAULT



# =====================================================================
# AGENT
# =====================================================================

def learning_strategy_choices():
    return [(x, x) for x in ReinventEnvironmentHelper.get_learning_strategies()]


class ReinventAgentTraining(TrainingStrategy):
    batch_size = models.IntegerField(default=64)
    unique_sequences = models.BooleanField(default=True)
    randomize_smiles = models.BooleanField(default=True)
    tb_isim = models.BooleanField(default=False)

    use_checkpoint = models.BooleanField(default=False)
    purge_memories = models.BooleanField(default=False)

    summary_csv_prefix = models.CharField(max_length=128, default="reinvent")

    learning_type = models.CharField(max_length=32, choices=learning_strategy_choices, default="dap")
    sigma = models.FloatField(default=128.0)
    rate = models.FloatField(default=0.0001)


class ReinventAgentValidation(ValidationStrategy):
    validate_every = models.IntegerField(default=50)
    validation_dataset = models.ForeignKey(ModelFile, null=True, blank=True, on_delete=models.SET_NULL)


class ReinventAgent(Model):
    environment = models.ForeignKey(ReinventEnvironment, on_delete=models.CASCADE)
    training = models.ForeignKey(ReinventAgentTraining, on_delete=models.PROTECT)
    validation = models.ForeignKey(ReinventAgentValidation, null=True, blank=True, on_delete=models.SET_NULL)
    output_model = models.ForeignKey(ModelFile, null=True, blank=True, on_delete=models.SET_NULL)

    # moved from Generator (keep defaults to avoid breaking existing configs)
    tb_logdir = models.CharField(max_length=255, default="tb_logs")
    json_out_config = models.CharField(max_length=255, default="_staged_learning.json")

    def getGenerator(self):
        # Mirrors DrugExAgent.getGenerator()
        return self.generator.order_by("-id").first()

    def _run_tag(self, gen) -> str:
        project_id = getattr(gen, "project_id", None)
        return f"project{project_id}_run{gen.pk}" if project_id else f"run{gen.pk}"

    # ------------------------------------------------------------------
    # Simple file helpers for staged learning
    # ------------------------------------------------------------------
    def _sl_dir(self) -> str:
        base = getattr(settings, "MEDIA_ROOT", ".")
        return os.path.join(base, "reinvent_sl")

    def get_toml_path(self, generator=None) -> str:
        gen = generator or self.getGenerator()
        if not gen:
            raise RuntimeError("This ReinventAgent is not attached to a Reinvent generator.")
        run_tag = self._run_tag(gen)
        return os.path.join(self._sl_dir(), f"staged_learning_{run_tag}.toml")

    def get_rl_log_path(self, generator=None) -> str:
        gen = generator or self.getGenerator()
        if not gen:
            raise RuntimeError("This ReinventAgent is not attached to a Reinvent generator.")
        run_tag = self._run_tag(gen)
        return os.path.join(self._sl_dir(), f"reinvent_rl_{run_tag}.log")

    def _results_prefix(self, gen) -> str:
        base = getattr(self.training, "summary_csv_prefix", None) or "reinvent"
        run_tag = self._run_tag(gen)
        return f"{base}_{run_tag}" if run_tag else base

    def get_csv_path(self, generator=None) -> str:
        """Resolve the staged-learning CSV path for this run.

        REINVENT writes one CSV per stage: {prefix}_1.csv, {prefix}_2.csv, …
        We want the **latest** stage (highest suffix number) because it
        contains the final RL results the user cares about.
        """
        gen = generator or self.getGenerator()
        if not gen:
            raise RuntimeError("This ReinventAgent is not attached to a Reinvent generator.")
        sl_dir = self._sl_dir()
        prefix = self._results_prefix(gen)
        legacy_prefix = getattr(self.training, "summary_csv_prefix", None) or "reinvent"

        import glob

        def _latest_numbered(pfx):
            """Return the CSV with the highest stage-number suffix, or None."""
            pattern = os.path.join(sl_dir, f"{pfx}_*.csv")
            matches = glob.glob(pattern)
            if not matches:
                return None
            # Extract the numeric suffix and pick the highest
            def _stage_num(p):
                base = os.path.basename(p)          # e.g. reinvent_project5_run43_2.csv
                stem = base.rsplit(".", 1)[0]        # reinvent_project5_run43_2
                parts = stem.rsplit("_", 1)          # ['reinvent_project5_run43', '2']
                try:
                    return int(parts[-1])
                except (ValueError, IndexError):
                    return 0
            matches.sort(key=_stage_num, reverse=True)
            return matches[0]

        # 1) Try numbered CSVs with the run-specific prefix (latest stage)
        latest = _latest_numbered(prefix)
        if latest:
            return latest

        # 2) Plain (un-numbered) file
        plain = os.path.join(sl_dir, f"{prefix}.csv")
        if os.path.isfile(plain):
            return plain

        # 3) Fallback: legacy prefix
        latest_legacy = _latest_numbered(legacy_prefix)
        if latest_legacy:
            return latest_legacy
        plain_legacy = os.path.join(sl_dir, f"{legacy_prefix}.csv")
        if os.path.isfile(plain_legacy):
            return plain_legacy

        # 4) Last resort: newest CSV matching prefix by mtime
        try:
            matches = glob.glob(os.path.join(sl_dir, f"{prefix}*.csv"))
            if not matches:
                matches = glob.glob(os.path.join(sl_dir, f"{legacy_prefix}*.csv"))
            if matches:
                matches.sort(key=lambda p: os.path.getmtime(p), reverse=True)
                return matches[0]
        except Exception:
            pass

        return os.path.join(sl_dir, f"{prefix}.csv")

    # ------------------------------------------------------------------
    # Build staged-learning TOML
    # ------------------------------------------------------------------
    def build_staged_toml(self, device="cuda:0", generator=None) -> str:
        gen = generator or self.getGenerator()
        if not gen:
            raise RuntimeError("This ReinventAgent is not attached to a Reinvent generator.")

        env = self.environment
        train_cfg = self.training

        # Generate unique tb_logdir with project and run IDs
        run_tag = self._run_tag(gen)
        tb_dir_path = os.path.join(self._sl_dir(), f"tb_logs_{run_tag}")

        lines = []
        lines.append('run_type = "staged_learning"')
        lines.append(f'device = "{device}"')
        lines.append(f'tb_logdir = "{tb_dir_path}"')
        lines.append(f'json_out_config = "{self.json_out_config}"')
        lines.append("")

        # PARAMETERS
        lines.append("[parameters]")
        lines.append(f"use_checkpoint = {str(train_cfg.use_checkpoint).lower()}")
        lines.append(f'summary_csv_prefix = "{self._results_prefix(gen)}"')

        # NOTE: model selection now supports Model OR ModelFile
        lines.append(f'prior_file = "{env.get_prior_path()}"')
        lines.append(f'agent_file = "{env.get_agent_path()}"')

        lines.append(f"batch_size = {train_cfg.batch_size}")
        lines.append(f"unique_sequences = {str(train_cfg.unique_sequences).lower()}")
        lines.append(f"randomize_smiles = {str(train_cfg.randomize_smiles).lower()}")
        lines.append(f"tb_isim = {str(train_cfg.tb_isim).lower()}")
        lines.append("")

        # LEARNING STRATEGY
        lines.append("[learning_strategy]")
        lines.append(f'type = "{train_cfg.learning_type}"')
        lines.append(f"sigma = {train_cfg.sigma}")
        lines.append(f"rate = {train_cfg.rate}")
        lines.append("")

        # DIVERSITY FILTER
        if env.diversity_filter:
            d = env.diversity_filter.to_reinvent()
            lines.append("[diversity_filter]")
            for k, v in d.items():
                lines.append(f'{k} = "{v}"' if isinstance(v, str) else f"{k} = {v}")
            lines.append("")

        # INCEPTION
        if env.inception_smiles:
            lines.append("[inception]")
            lines.append(f'smiles_file = "{env.inception_smiles.file.path}"')
            if env.inception_memory_size:
                lines.append(f"memory_size = {env.inception_memory_size}")
            if env.inception_sample_size:
                lines.append(f"sample_size = {env.inception_sample_size}")
            lines.append("")

        # STAGES
        stages = list(gen.stages.order_by("order"))
        if not stages:
            raise RuntimeError(
                "This run has no stages defined. "
                "Go to the Multi-Stage Learning tab, select this run, "
                "assign scoring components, and create at least one stage."
            )

        sl_dir = self._sl_dir()
        os.makedirs(sl_dir, exist_ok=True)
        run_tag = self._run_tag(gen)

        for stage_idx, st in enumerate(stages):
            lines.append("[[stage]]")

            default_chkpt = os.path.join(sl_dir, f"agent_stage{st.order}_{run_tag}.chkpt")
            chkpt_path = st.resolve_checkpoint_path(default_chkpt)
            lines.append(f'chkpt_file = "{chkpt_path}"')

            lines.append(f'termination = "{st.termination_type}"')
            lines.append(f"max_score = {st.max_score}")
            lines.append(f"min_steps = {st.min_steps}")
            lines.append(f"max_steps = {st.max_steps}")
            lines.append("")

            prefix = "stage.scoring"

            # EXTERNAL SCORING FILE
            if st.scoring_source == "file":
                if not st.scoring_file:
                    raise RuntimeError("Stage scoring_source='file' but scoring_file is not set.")
                agg = st.aggregation_type or "geometric_mean"
                lines.append(f"[{prefix}]")
                lines.append(f'type = "{agg}"')
                lines.append(f'filename = "{st.scoring_file.file.path}"')
                lines.append('filetype = "toml"')
                lines.append("")
                continue

            # INLINE SCORING — use stage aggregation_type (defaults to geometric_mean)
            agg = st.aggregation_type or "geometric_mean"

            # Every stage automatically gets the default custom_alerts
            # component, so there is always at least one scoring component.

            lines.append(f"[{prefix}]")
            lines.append(f'type = "{agg}"')
            lines.append("")

            # Translate any legacy / RDKit-style property names stored in the DB
            # to the canonical REINVENT 4 component names from Table 2.
            # Keys  = names that may exist in old PropertyScorer.property_name rows.
            # Values = correct component names recognised by REINVENT 4 (case-
            #          insensitive lookup, but we use canonical capitalisation).
            _PROP_NAME_MAP = {
                # ── QED ─────────────────────────────────────────────────
                "QED": "Qed",
                "qed": "Qed",
                # ── Lipophilicity ────────────────────────────────────────
                "MolLogP": "SlogP",
                "logP": "SlogP",
                "LogP": "SlogP",
                # ── Molecular weight ─────────────────────────────────────
                "MolWt": "MolecularWeight",
                "MW": "MolecularWeight",
                "mol_weight": "MolecularWeight",
                # ── HBond donors / acceptors ─────────────────────────────
                "NumHDonors": "HBondDonors",
                "NumHAcceptors": "HBondAcceptors",
                # ── Rotatable bonds ──────────────────────────────────────
                "NumRotatableBonds": "NumRotBond",
                # ── Rings — old RDKit names that have no direct equivalent
                #    are mapped to the closest valid REINVENT 4 component ──
                "NumSaturatedRings":       "NumRings",
                "NumAromaticHeterocycles": "NumAromaticRings",
                "NumAliphaticHeterocycles":"NumAliphaticRings",
                "NumSaturatedHeterocycles":"NumRings",
                "NumAromaticCarbocycles":  "NumAromaticRings",
                "NumAliphaticCarbocycles": "NumAliphaticRings",
                "NumSaturatedCarbocycles": "NumRings",
                # ── sp hybridisation ────────────────────────────────────
                # REINVENT 4 registers these lowercase (numsp, numsp2, numsp3)
                "Numsp":  "numsp",
                "Numsp2": "numsp2",
                "Numsp3": "numsp3",
                # ── Stereocenters ────────────────────────────────────────
                "NumStereocenters":    "NumAtomStereoCenters",
                "NumStereoCenters":    "NumAtomStereoCenters",
                "num_stereocenters":   "NumAtomStereoCenters",
                # ── Graph / topology ─────────────────────────────────────
                "graph_length": "GraphLength",
                # ── sp3 fraction ─────────────────────────────────────────
                "FractionCSP3": "Csp3",
                "FCsp3":        "Csp3",
                # ── New atom/ring counts ──────────────────────────────────────
                "HeavyAtomCount":       "NumHeavyAtoms",
                "num_heavy_atoms":      "NumHeavyAtoms",
                "NumHeteroatoms":       "NumHeteroAtoms",
                "num_heteroatoms":      "NumHeteroAtoms",
                "NumAliphRings":        "NumAliphaticRings",
                "num_aliphatic_rings":  "NumAliphaticRings",
                "NumAtomStereoCenters": "NumAtomSteroCenters",
                "num_stereo_centers":   "NumAtomSteroCenters",
            }

            # Default component_params for REINVENT 4 components that
            # require endpoint-level params.  When the user has not supplied
            # explicit values we fall back to these sensible defaults so the
            # TOML is always valid.
            # IMPORTANT: REINVENT's collect_params() aggregates per-
            # endpoint param values into lists. That means every value
            # written here must be a **scalar** (string / int / bool)
            # because collect_params will wrap it in a list.
            #
            # The ONE exception is TanimotoSimilarity/Distance where the
            # Pydantic field is List[List[str]] for smiles, so the
            # endpoint-level value must already be a list of strings
            # (collect_params wraps it into the outer list).
            _COMPONENT_DEFAULT_PARAMS = {
                "TanimotoSimilarity": {
                    "smiles": ["c1ccccc1"],   # List[str] → collect_params → List[List[str]]
                    "radius": 3,
                    "use_counts": True,
                    "use_features": False,
                },
                "TanimotoDistance": {
                    "smiles": ["c1ccccc1"],   # List[str] → collect_params → List[List[str]]
                    "radius": 3,
                    "use_counts": True,
                    "use_features": False,
                },
                "GroupCount": {
                    "smarts": "[#6]",         # scalar → collect_params → List[str]
                },
                "MatchingSubstructure": {
                    "smarts": "[#6]",         # scalar → collect_params → List[str]
                    "use_chirality": False,
                },
                "PMI": {
                    "property": "npr1",       # scalar → collect_params → List[str]
                },
            }

            # Property scorers attached to this stage
            for ps in st.property_scorers.order_by("id"):
                prop = _PROP_NAME_MAP.get(ps.property_name, ps.property_name)
                lines.append(f"[[{prefix}.component]]")
                lines.append(f"[{prefix}.component.{prop}]")
                lines.append(f"[[{prefix}.component.{prop}.endpoint]]")
                lines.append(f'name = "{ps.name}"')
                lines.append(f"weight = {ps.weight}")

                # Emit component_params as endpoint-level params.
                # Merge user-supplied params over defaults so required
                # fields are always present.
                defaults = _COMPONENT_DEFAULT_PARAMS.get(prop, {})
                cparams = dict(defaults)  # start with defaults
                if ps.component_params and isinstance(ps.component_params, dict):
                    cparams.update(ps.component_params)  # user overrides

                if cparams:
                    # Keys whose endpoint-level value must remain a list
                    # (because the Pydantic field is List[List[...]]).
                    _LIST_KEYS = {"smiles"}

                    for pk, pv in cparams.items():
                        # Unwrap single-element lists for params that
                        # should be scalar.  collect_params() will wrap
                        # them back into a list later.
                        if (isinstance(pv, list) and len(pv) == 1
                                and pk not in _LIST_KEYS):
                            pv = pv[0]

                        if isinstance(pv, list):
                            items = ", ".join(f'"{x}"' if isinstance(x, str) else str(x) for x in pv)
                            lines.append(f"params.{pk} = [{items}]")
                        elif isinstance(pv, bool):
                            lines.append(f"params.{pk} = {str(pv).lower()}")
                        elif isinstance(pv, str):
                            lines.append(f'params.{pk} = "{pv}"')
                        else:
                            lines.append(f"params.{pk} = {pv}")

                tr = ps.build_transform()
                if tr:
                    lines.append(f"[{prefix}.component.{prop}.endpoint.transform]")
                    for k, v in tr.items():
                        lines.append(f'{k} = "{v}"' if isinstance(v, str) else f"{k} = {v}")
                lines.append("")

            # QSAR model scorers attached to this stage
            for ms in st.model_scorers.order_by("id"):
                lines.append(f"[[{prefix}.component]]")
                lines.append(f"[{prefix}.component.qsar_models]")
                lines.append(f"[[{prefix}.component.qsar_models.endpoint]]")
                lines.append(f'name = "{ms.name}"')
                lines.append(f"weight = {ms.weight}")
                tr = ms.build_transform()
                if tr:
                    lines.append(f"[{prefix}.component.qsar_models.endpoint.transform]")
                    for k, v in tr.items():
                        lines.append(f'{k} = "{v}"' if isinstance(v, str) else f"{k} = {v}")
                lines.append("")

            # SMARTS scorers attached to this stage (user-defined UnwantedSmartsScorer)
            has_user_smarts = False
            for us in st.smarts_scorers.order_by("id"):
                if not us.enabled:
                    continue
                has_user_smarts = True
                patterns = us.load_patterns()

                lines.append(f"[[{prefix}.component]]")
                lines.append(f"[{prefix}.component.custom_alerts]")
                lines.append(f"[[{prefix}.component.custom_alerts.endpoint]]")
                lines.append(f'name = "{us.name}"')
                lines.append(f"weight = {us.weight}")
                lines.append("params.smarts = [")
                for p in patterns:
                    lines.append(f'  "{p}",')
                lines.append("]")
                lines.append("")


            # ── Default bad-SMARTS penalty (first stage only, no user SMARTS) ─
            # REINVENT applies custom_alerts once per run; including it
            # only in the first stage avoids redundant duplication.
            if stage_idx == 0 and not has_user_smarts:
                _bad_smarts_w = getattr(gen, "bad_smarts_weight", 1.0)
                if _bad_smarts_w is None:
                    _bad_smarts_w = 1.0
                lines.append(f"[[{prefix}.component]]")
                lines.append(f"[{prefix}.component.custom_alerts]")
                lines.append(f"[[{prefix}.component.custom_alerts.endpoint]]")
                lines.append('name = "Unwanted SMARTS (default)"')
                lines.append(f"weight = {_bad_smarts_w}")
                lines.append("params.smarts = [")
                for pat in UNWANTED_SMARTS_DEFAULT:
                    lines.append(f'  "{pat}",')
                lines.append("]")
                lines.append("")

            lines.append("")

        body = "\n".join(lines) + "\n"
        toml_path = self.get_toml_path(generator=gen)
        with open(toml_path, "w", encoding="utf-8") as fh:
            fh.write(body)

        return toml_path

    # ------------------------------------------------------------------
    # RUN STAGED LEARNING
    # ------------------------------------------------------------------
    def run_staged_learning(self, device="cuda:0", generator=None) -> str:
        gen = generator or self.getGenerator()
        if not gen:
            raise RuntimeError("This ReinventAgent is not attached to a Reinvent generator.")

        run_tag = self._run_tag(gen)
        print(f"[ReinventAgent.run_staged_learning] Generator={gen.id}, ProjectID={gen.project_id}, RunTag={run_tag}")

        toml_path = self.build_staged_toml(device=device, generator=gen)

        reinvent_bin = (
            getattr(settings, "REINVENT_BIN", None)
            or os.environ.get("REINVENT_BIN")
            or shutil.which("reinvent")
        )
        if not reinvent_bin:
            raise RuntimeError("REINVENT binary not found. Set REINVENT_BIN in settings or env.")

        cmd = [reinvent_bin, toml_path]
        rl_dir = self._sl_dir()
        os.makedirs(rl_dir, exist_ok=True)

        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            cwd=rl_dir,
        )

        lines = [ln for ln in (proc.stdout or [])]
        rc = proc.wait()

        log_path = self.get_rl_log_path(generator=gen)
        log_text = f"[CMD] {' '.join(cmd)}\n{''.join(lines)}"
        with open(log_path, "w", encoding="utf-8") as fh:
            fh.write(log_text)

        if rc != 0:
            raise RuntimeError(f"REINVENT staged learning failed (exit {rc}). See log: {log_path}")

        # Persist results CSV into the run's ModelFile (if available)
        try:
            csv_path = self.get_csv_path(generator=gen)
            if os.path.isfile(csv_path):
                with open(csv_path, "rb") as fh:
                    data = fh.read()
                result_file = gen.get_results_file()
                _overwrite_filefield(result_file, data, filename=os.path.basename(result_file.file.name))
        except Exception:
            pass

        return toml_path


class Reinvent(Generator):
    # keep fields if you already migrated data; frontend can keep using them
    environment = models.ForeignKey(ReinventEnvironment, on_delete=models.CASCADE)
    agent = models.ForeignKey(ReinventAgent, on_delete=models.PROTECT, related_name="generator")
    bad_smarts_weight = models.FloatField(default=1.0, help_text="Weight for the default bad-SMARTS penalty (0–1)")

    RESULTS_FILE_NOTE = "reinvent_sl_results"

    def _run_tag(self) -> str:
        return f"project{self.project_id}_run{self.pk}" if self.project_id else f"run{self.pk}"

    def get_results_file(self) -> ModelFile:
        """Return a ModelFile for staged-learning CSV results."""
        mf = self.files.filter(kind=ModelFile.AUXILIARY, note=self.RESULTS_FILE_NOTE).first()
        if mf is None:
            filename = f"reinvent_{self._run_tag()}.csv"
            mf = ModelFile.create(
                self,
                filename,
                ContentFile(b""),
                note=self.RESULTS_FILE_NOTE,
            )
        return mf

    def _extract_smiles_from_csv(self, csv_path: str) -> list[str]:
        """Extract SMILES from CSV file.

        Returns list of SMILES strings.
        """
        smiles = []
        with open(csv_path, "r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh)
            fieldnames = reader.fieldnames or []
            smiles_key = None
            for name in fieldnames:
                if name and name.strip().lower() == "smiles":
                    smiles_key = name
                    break
            if not smiles_key:
                return smiles
            for row in reader:
                val = (row.get(smiles_key) or "").strip()
                if val:
                    smiles.append(val)
        return smiles

    def _extract_smiles_with_scores_from_csv(self, csv_path: str, min_score: float = None) -> list[tuple[str, float]]:
        """Extract SMILES and scores from CSV file with optional filtering.

        Returns list of tuples (smiles, score).

        This method supports multiple CSV formats:
        1. REINVENT 4.x format with headers: Agent,Prior,Target,Score,SMILES,SMILES_state,...
        2. Alternative REINVENT format: step,agent,smiles,score,component1,...
        3. Legacy format without headers: component1_score,component2_score,...,total_score,SMILES
        """
        import logging
        logger = logging.getLogger(__name__)

        results = []
        with open(csv_path, "r", encoding="utf-8", newline="") as fh:
            reader = csv.reader(fh)
            row_number = 0
            header_row = None
            smiles_col_idx = None
            score_col_idx = None

            for row in reader:
                row_number += 1
                if len(row) < 2:
                    continue

                # Detect header row
                if row_number == 1:
                    # Check if first row is a header
                    if any(h.lower() in ['smiles', 'step', 'agent', 'score'] for h in row):
                        header_row = [h.lower() for h in row]
                        # Find SMILES and score columns
                        # Note: Look for exact column name match first, then check for substring
                        for idx, h in enumerate(header_row):
                            h_stripped = h.strip()
                            # Find SMILES column - prefer exact 'smiles' match over 'smiles_state'
                            if h_stripped == 'smiles':
                                smiles_col_idx = idx
                            # Find score column - prefer exact 'score' match
                            if h_stripped in ['score', 'total_score', 'avg_score']:
                                score_col_idx = idx
                        logger.info(f"[CSV Parse] Detected header row with SMILES column at index {smiles_col_idx}, score at {score_col_idx}")
                        continue

                try:
                    smiles = None
                    total_score = None

                    # Extract SMILES and score based on detected format
                    if smiles_col_idx is not None and score_col_idx is not None:
                        # Use detected column indices
                        if len(row) <= max(smiles_col_idx, score_col_idx):
                            continue
                        smiles = row[smiles_col_idx].strip()
                        total_score = float(row[score_col_idx])
                    else:
                        # No header detected - try to auto-detect SMILES column
                        # SMILES should be a string containing chemical notation
                        # Scan columns from right to left to find the SMILES column
                        for idx in range(len(row) - 1, -1, -1):
                            val = row[idx].strip()
                            if not val:
                                continue

                            # Try to identify SMILES by checking for typical characters
                            if any(c in val for c in ['C', 'N', 'O', 'S', 'c', 'n', 'o', 's', '(', ')', '=', '#', '[', ']']):
                                # Check if it's not just a number
                                try:
                                    float(val)
                                    continue  # It's a number, not SMILES
                                except ValueError:
                                    # Good, it's not a plain number
                                    smiles = val
                                    # Score should be in a column before SMILES
                                    for score_idx in range(idx - 1, -1, -1):
                                        try:
                                            total_score = float(row[score_idx])
                                            break
                                        except (ValueError, IndexError):
                                            continue
                                    break

                        # If still not found, fall back to legacy format: score in second-to-last, SMILES in last
                        if smiles is None and len(row) >= 2:
                            try:
                                total_score = float(row[-2])
                                smiles = row[-1].strip()
                            except (ValueError, IndexError):
                                continue

                    if not smiles or total_score is None:
                        continue

                    # Validate that SMILES is not just a number (likely an index)
                    try:
                        float(smiles)
                        # If we can convert SMILES to float, it's probably a row index, skip it
                        if row_number <= 10:  # Only log first 10 warnings
                            logger.warning(f"[CSV Parse] Row {row_number}: Skipping potential index '{smiles}' instead of SMILES")
                        continue
                    except ValueError:
                        # Good, it's not a plain number, proceed
                        pass

                    # Basic SMILES validation - should contain typical SMILES characters
                    if not any(c in smiles for c in ['C', 'N', 'O', 'S', 'c', 'n', 'o', 's', '(', ')', '=', '#']):
                        if row_number <= 10:  # Only log first 10 warnings
                            logger.warning(f"[CSV Parse] Row {row_number}: Skipping invalid SMILES '{smiles}'")
                        continue

                    # Filter by minimum score if specified
                    if min_score is not None and total_score < min_score:
                        continue

                    results.append((smiles, total_score))
                except (ValueError, IndexError) as e:
                    if row_number <= 10:  # Only log first 10 errors
                        logger.debug(f"[CSV Parse] Row {row_number}: Failed to parse - {e}")
                    continue

        return results


    def build_staged_toml(self, device="cuda:0") -> str:
        return self.agent.build_staged_toml(device=device, generator=self)

    def run_staged_learning(self, device="cuda:0") -> str:
        return self.agent.run_staged_learning(device=device, generator=self)

    def get(self, n_samples, min_score=None) -> list[str]:
        """
        Get n_samples SMILES strings from staged learning results.

        This method is used by the Compounds -> Generated Sets interface
        to retrieve molecules from the staged learning run results.

        For Reinvent staged learning, molecules are already generated during
        the run and stored in CSV files. This method reads from those CSV files
        instead of running new sampling.

        :param n_samples: Number of SMILES to return
        :param min_score: Optional minimum score threshold for filtering
        :return: List of SMILES strings
        """
        import logging

        logger = logging.getLogger(__name__)
        logger.info(f"[Reinvent.get] Getting {n_samples} SMILES from staged learning results (min_score={min_score})")

        # Get the CSV path for this run
        csv_path = self.agent.get_csv_path(generator=self)


        logger.info(f"[Reinvent.get] Reading SMILES from CSV: {csv_path}")

        # Check if file exists and log its size
        if os.path.exists(csv_path):
            file_size = os.path.getsize(csv_path)
            logger.info(f"[Reinvent.get] CSV file size: {file_size} bytes")

            # Log first few lines for debugging
            try:
                with open(csv_path, 'r') as f:
                    first_lines = []
                    for i in range(5):
                        line = f.readline()
                        if not line:
                            break
                        first_lines.append(line.strip())
                    logger.info(f"[Reinvent.get] First {len(first_lines)} lines of CSV:")
                    for i, line in enumerate(first_lines):
                        logger.info(f"[Reinvent.get]   Line {i+1}: {line[:200] if len(line) > 200 else line}")
            except Exception as e:
                logger.warning(f"[Reinvent.get] Could not read CSV preview: {e}")
        else:
            logger.error(f"[Reinvent.get] CSV file does not exist!")

        # Extract SMILES with scores from the CSV file
        smiles_with_scores = self._extract_smiles_with_scores_from_csv(csv_path, min_score=min_score)

        if not smiles_with_scores:
            # Get statistics about available scores to help user
            all_smiles_with_scores = self._extract_smiles_with_scores_from_csv(csv_path, min_score=None)

            if min_score is not None and all_smiles_with_scores:
                # Calculate score statistics
                all_scores = [score for _, score in all_smiles_with_scores]
                max_available = max(all_scores)
                min_available = min(all_scores)
                median_score = sorted(all_scores)[len(all_scores) // 2]

                logger.error(
                    f"[Reinvent.get] No molecules found with score >= {min_score}. "
                    f"Available score range: {min_available:.4f} to {max_available:.4f} "
                    f"(median: {median_score:.4f})"
                )

                raise ValueError(
                    f"No molecules found with score >= {min_score}. "
                    f"The highest score in this run is {max_available:.4f}. "
                    f"Available score range: {min_available:.4f} to {max_available:.4f} (median: {median_score:.4f}). "
                    f"Please lower the minimum score threshold to at least {max_available:.4f} or less."
                )
            else:
                # No SMILES found at all - CSV is malformed or empty
                logger.error(
                    f"[Reinvent.get] No valid SMILES found in CSV file. "
                    f"The file may be empty, corrupted, or in an unexpected format. "
                    f"Please check the CSV file at: {csv_path}"
                )
                raise ValueError(
                    f"No SMILES found in CSV file {csv_path}. "
                    f"The CSV file may be empty or improperly formatted. "
                    f"Expected format: CSV with 'smiles' and 'score' columns (REINVENT 4.x), "
                    f"or legacy format with score in second-to-last column and SMILES in last column."
                )

        logger.info(f"[Reinvent.get] Found {len(smiles_with_scores)} SMILES matching criteria (min_score={min_score})")

        # Check if we have fewer molecules than requested
        if len(smiles_with_scores) < n_samples:
            logger.warning(
                f"[Reinvent.get] Only {len(smiles_with_scores)} molecules available "
                f"(requested {n_samples}). Returning all available molecules."
            )
            if min_score is not None:
                logger.warning(
                    f"[Reinvent.get] Consider lowering the score threshold (current: {min_score}) "
                    f"to get more molecules."
                )

        # Sort by score descending to get best molecules first
        smiles_with_scores.sort(key=lambda x: x[1], reverse=True)

        # Log score range
        if smiles_with_scores:
            max_score = smiles_with_scores[0][1]
            min_found_score = smiles_with_scores[-1][1]
            logger.info(f"[Reinvent.get] Score range: {min_found_score:.4f} to {max_score:.4f}")

        # Extract just the SMILES strings
        smiles_list = [smi for smi, score in smiles_with_scores]

        # Return requested number of SMILES (or all if fewer available)
        actual_count = min(len(smiles_list), n_samples)
        result = smiles_list[:actual_count]
        logger.info(f"[Reinvent.get] Returning {len(result)} SMILES (requested {n_samples})")

        return result

    def __str__(self):
        """String representation for the model."""
        # Try different ways to get a meaningful name
        name = getattr(self, 'name', None)
        if name and name.strip():
            return name

        pk = getattr(self, 'pk', None)
        project_id = getattr(self, 'project_id', None)

        if pk:
            if project_id:
                return f"Reinvent (Project {project_id}, Run {pk})"
            return f"Reinvent Run {pk}"

        return "Reinvent (unsaved)"

# =====================================================================
# 6. STAGE MODEL
# =====================================================================

class ReinventStage(models.Model):
    generator = models.ForeignKey("Reinvent", related_name="stages", on_delete=models.CASCADE)
    order = models.IntegerField()

    # existing
    chkpt_file = models.ForeignKey(
        ModelFile, null=True, blank=True, on_delete=models.SET_NULL,
        related_name="reinvent_stage_checkpoints"
    )

    # NEW: pick an existing trained net as the checkpoint source
    chkpt_net = models.ForeignKey(
        "ReinventNet", null=True, blank=True, on_delete=models.SET_NULL,
        related_name="reinvent_stage_checkpoint_nets"
    )

    termination_type = models.CharField(max_length=64, default="simple")
    max_score = models.FloatField(default=1.0)
    min_steps = models.IntegerField(default=1)
    max_steps = models.IntegerField(default=100)

    scoring_source = models.CharField(
        max_length=16,
        choices=[("inline", "Inline"), ("file", "File")],
        default="inline"
    )

    aggregation_type = models.CharField(
        max_length=64,
        choices=[
            ("geometric_mean", "Geometric Mean (Balanced - all scores must be good)"),
            ("arithmetic_mean", "Arithmetic Mean (Flexible - based on weights)"),
        ],
        default="geometric_mean",
        null=True, blank=True
    )

    scoring_file = models.ForeignKey(
        ModelFile,
        null=True, blank=True,
        on_delete=models.SET_NULL,
        related_name="reinvent_stage_scoring_files"
    )

    # Scoring components assigned to this stage
    property_scorers = models.ManyToManyField(
        "PropertyScorer", blank=True, related_name="stages"
    )
    model_scorers = models.ManyToManyField(
        "GenUIModelScorer", blank=True, related_name="stages"
    )
    smarts_scorers = models.ManyToManyField(
        "UnwantedSmartsScorer", blank=True, related_name="stages"
    )

    class Meta:
        ordering = ["order"]

    def resolve_checkpoint_path(self, default_path: str) -> str:
        if self.chkpt_net:
            return self.chkpt_net.get_active_checkpoint_path()
        if self.chkpt_file:
            return self.chkpt_file.file.path
        return default_path

# =====================================================================
# 7. PERFORMANCE LOGGING
# =====================================================================



class ModelPerformanceReinvent(ModelPerformance):
    epoch = models.IntegerField()
    step = models.IntegerField()
    created = models.DateTimeField(default=timezone.now, db_index=True)

    isOnValidationSet = models.BooleanField(default=False)
    note = models.CharField(max_length=128, blank=True)

    stage_index = models.IntegerField(null=True, blank=True)
    avg_score = models.FloatField(null=True, blank=True)
    fraction_valid = models.FloatField(null=True, blank=True)
    avg_nll = models.FloatField(null=True, blank=True)
    unique_scaffolds = models.IntegerField(null=True, blank=True)

