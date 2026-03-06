# genui/generators/extensions/genuireinvent/views.py

from __future__ import annotations

import logging
import os

from django.conf import settings
from django.db import close_old_connections

from rest_framework import status, viewsets
from rest_framework.decorators import action
from rest_framework.response import Response

from genui.models.views import ModelViewSet, MetricsViewSet
from genui.models.models import ModelFile, AlgorithmMode
from genui.models.serializers import ModelFileSerializer

from . import models, serializers
from .genuimodels import builders
from .tasks import buildReinventModel, runReinventStagedLearning
from .genuisetup import ensure_reinvent_prior
from genui.utils.extensions.tasks.utils import runTask

log = logging.getLogger(__name__)


def _auto_detect_device(requested: str = "auto") -> str:
    """Resolve device string. If 'auto', pick cuda:0 when available, else cpu."""
    if requested and requested != "auto":
        return requested
    try:
        import torch
        if torch.cuda.is_available():
            return "cuda:0"
    except ImportError:
        pass
    return "cpu"


class ReinventNetViewSet(ModelViewSet):
    """
    CRUD for ReinventNet + action to run REINVENT's preprocessor.
    Artifacts are stored as AUX ModelFiles (hashed paths under media/).
    """
    queryset = models.ReinventNet.objects.order_by("-created")
    serializer_class = serializers.ReinventNetSerializer
    init_serializer_class = serializers.ReinventNetInitSerializer
    owner_relation = "project__owner"
    builder_class = builders.ReinventNetBuilder
    build_task = buildReinventModel

    def get_builder_kwargs(self):
        return {"model_class": models.ReinventNet.__name__}

    @action(detail=True, methods=["post"], url_path="prepare-corpus")
    def prepare_corpus(self, request, pk=None):
        """
        POST /reinvent/networks/{id}/prepare-corpus/
        Runs datapipeline, writes cleaned corpus + train/valid splits to AUX files,
        and returns counts + file references.
        """
        try:
            net = self.get_object()
            train_mf, valid_mf = net.prepareData()

            def _count_lines(path: str) -> int:
                try:
                    with open(path, "r", encoding="utf-8") as fh:
                        return sum(1 for ln in fh if ln.strip())
                except FileNotFoundError:
                    return 0

            return Response(
                {
                    "prepared_train": _count_lines(train_mf.path),
                    "prepared_valid": _count_lines(valid_mf.path),
                    "train_file": train_mf.path,
                    "valid_file": valid_mf.path,
                    "preview_file": net.corpusPreviewFile.path,
                    "full_file": net.corpusFullFile.path,
                    "split_method": getattr(getattr(net, "validationStrategy", None), "split_method", None),
                },
                status=status.HTTP_201_CREATED,
            )
        except Exception as e:
            log.exception("prepare_corpus failed for ReinventNet pk=%s", pk)
            return Response({"error": repr(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

    @action(detail=False, methods=["get"], url_path="prior-status")
    def prior_status(self, request):
        try:
            prior_path = models._resolve_reinvent_prior_path()
            return Response(
                {
                    "resolved_path": prior_path,
                    "exists": bool(prior_path and os.path.isfile(prior_path)),
                    "error": None,
                }
            )
        except Exception as e:
            return Response(
                {
                    "resolved_path": None,
                    "exists": False,
                    "error": repr(e),
                }
            )

    @action(detail=False, methods=["post"], url_path="resolve-prior")
    def resolve_prior(self, request):
        """
        POST /reinvent/networks/resolve-prior/
        Downloads and verifies the REINVENT prior file.
        """
        try:
            force = request.data.get("force", False) if request.data else False
            prior_path = ensure_reinvent_prior(force=force)
            return Response(
                {
                    "status": "success",
                    "resolved_path": prior_path,
                    "exists": os.path.isfile(prior_path),
                    "message": f"REINVENT prior resolved successfully at {prior_path}",
                },
                status=status.HTTP_200_OK,
            )
        except Exception as e:
            log.exception("resolve_prior failed")
            return Response(
                {
                    "status": "error",
                    "resolved_path": None,
                    "exists": False,
                    "message": str(e),
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

    @action(detail=True, methods=["get"], url_path="training-log")
    def training_log(self, request, pk=None):
        """
        GET /reinvent/networks/{id}/training-log/
        Returns the training log content for this ReinventNet.
        """
        try:
            net = self.get_object()
            log_file = net.trainLogFile

            if not log_file.file:
                return Response(
                    {"content": "", "exists": False, "message": "Training log not available yet"},
                    status=status.HTTP_200_OK,
                )

            try:
                with log_file.file.open("r") as f:
                    content = f.read()
                return Response(
                    {"content": content, "exists": True, "path": log_file.file.name},
                    status=status.HTTP_200_OK,
                )
            except Exception as read_err:
                log.warning(f"Could not read training log: {read_err}")
                return Response(
                    {"content": "", "exists": False, "message": f"Error reading log: {str(read_err)}"},
                    status=status.HTTP_200_OK,
                )
        except Exception as e:
            log.exception("training_log failed for ReinventNet pk=%s", pk)
            return Response({"error": repr(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

    @action(detail=True, methods=["get"], url_path="training-tb")
    def training_tb(self, request, pk=None):
        """
        GET /reinvent/networks/{id}/training-tb/
        Returns summary statistics from TensorBoard logs for this ReinventNet.
        """
        try:
            net = self.get_object()
            tb_dir = os.path.join(settings.MEDIA_ROOT, "models", f"tb_TL_{net.pk}")
            if not os.path.isdir(tb_dir):
                return Response(
                    {"exists": False, "message": "TensorBoard log dir not found", "tb_dir": tb_dir},
                    status=status.HTTP_200_OK,
                )

            try:
                from tensorboard.backend.event_processing.event_accumulator import EventAccumulator
            except Exception as import_err:
                return Response(
                    {"exists": False, "message": f"TensorBoard reader unavailable: {str(import_err)}"},
                    status=status.HTTP_200_OK,
                )

            def _read_scalars(path: str):
                try:
                    ea = EventAccumulator(path)
                    ea.Reload()
                    tags = ea.Tags().get("scalars", [])
                    return ea, tags
                except Exception:
                    return None, []

            def _get_series(ea, tags, tag_candidates):
                if ea is None:
                    return None, []
                for tag in tag_candidates:
                    if tag in tags:
                        vals = ea.Scalars(tag)
                        if vals:
                            return tag, vals
                return None, []

            def _latest(vals):
                if not vals:
                    return None
                last = vals[-1]
                return {"step": last.step, "value": last.value, "wall_time": last.wall_time}

            def _best(vals):
                if not vals:
                    return None
                best = min(vals, key=lambda x: x.value)
                return {"step": best.step, "value": best.value, "wall_time": best.wall_time}

            def _series(vals, limit=500):
                if not vals:
                    return []
                tail = vals[-limit:] if limit and len(vals) > limit else vals
                return [{"step": v.step, "value": v.value, "wall_time": v.wall_time} for v in tail]

            # Root TB scalars
            ea_root, tags_root = _read_scalars(tb_dir)
            valid_tag, valid_vals = _get_series(ea_root, tags_root, ["valid/nll", "validation/nll", "valid/loss", "validation/loss"])
            train_tag, train_vals = _get_series(ea_root, tags_root, ["train/nll", "training/nll", "train/loss", "training/loss"])

            # Subfolder TB scalars (Sample/Training/Validation Loss)
            subfolders = {
                "sample": "A_Mean NLL loss_Sample Loss",
                "train": "A_Mean NLL loss_Training Loss",
                "valid": "A_Mean NLL loss_Validation Loss",
            }
            sub_scalars = {}
            for key, sub in subfolders.items():
                sub_path = os.path.join(tb_dir, sub)
                ea_sub, tags_sub = _read_scalars(sub_path) if os.path.isdir(sub_path) else (None, [])
                tag, vals = _get_series(ea_sub, tags_sub, tags_sub)
                sub_scalars[key] = {
                    "dir": sub_path,
                    "tag": tag,
                    "latest": _latest(vals),
                    "best": _best(vals),
                    "series": _series(vals),
                }

            payload = {
                "exists": True,
                "tb_dir": tb_dir,
                "tags": tags_root,
                "valid": {
                    "tag": valid_tag,
                    "latest": _latest(valid_vals),
                    "best": _best(valid_vals),
                    "series": _series(valid_vals),
                },
                "train": {
                    "tag": train_tag,
                    "latest": _latest(train_vals),
                    "best": _best(train_vals),
                    "series": _series(train_vals),
                },
                "subfolders": sub_scalars,
            }

            return Response(payload, status=status.HTTP_200_OK)
        except Exception as e:
            log.exception("training_tb failed for ReinventNet pk=%s", pk)
            return Response({"error": repr(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)


class ReinventDiversityFilterViewSet(viewsets.ModelViewSet):
    queryset = models.ReinventDiversityFilter.objects.all().order_by("-id")
    serializer_class = serializers.ReinventDiversityFilterSerializer


    @action(detail=False, methods=["get"], url_path="available-types")
    def available_types(self, request):
        try:
            types_ = models.ReinventEnvironmentHelper.get_diversity_filters()
            return Response({"types": types_})
        except Exception as e:
            log.exception("available_types failed")
            return Response({"error": repr(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)


class ScoreModifierViewSet(viewsets.ModelViewSet):
    queryset = models.ScoreModifier.objects.all()
    serializer_class = serializers.ScoreModifierSerializer

    @action(detail=False, methods=["post"], url_path="test")
    def test(self, request):
        """
        POST /score-modifiers/test/
        Test a score modifier with example inputs to preview its transformation.

        Request body:
        {
            "inputs": [0.1, 0.2, 0.3, ...],
            "params": {"type": "ClippedScore", "upper": 1.0, "lower": 0.0, ...}
        }

        Returns:
        {
            "results": [modified_value1, modified_value2, ...]
        }
        """
        try:
            inputs = request.data.get('inputs', [])
            params = request.data.get('params', {})

            if not inputs:
                return Response(
                    {"error": "No inputs provided"},
                    status=status.HTTP_400_BAD_REQUEST
                )

            modifier_type = params.get('type', 'ClippedScore')

            # Apply modifier transformation
            results = []
            for input_value in inputs:
                try:
                    value = float(input_value)

                    if modifier_type == 'ClippedScore':
                        # Clipped score transformation
                        upper = float(params.get('upper', 1.0))
                        lower = float(params.get('lower', 0.0))
                        high = float(params.get('high', 1.0))
                        low = float(params.get('low', 0.0))
                        smooth = params.get('smooth', True)

                        if value >= upper:
                            result = high
                        elif value <= lower:
                            result = low
                        elif smooth:
                            # Smooth sigmoid transformation
                            import math
                            # Normalize to [0, 1]
                            normalized = (value - lower) / (upper - lower) if upper != lower else 0.5
                            # Apply sigmoid
                            sigmoid = 1 / (1 + math.exp(-10 * (normalized - 0.5)))
                            # Scale to [low, high]
                            result = low + (high - low) * sigmoid
                        else:
                            # Linear interpolation
                            result = low + (high - low) * ((value - lower) / (upper - lower))

                    elif modifier_type == 'SmoothHump':
                        # Smooth hump (Gaussian) transformation
                        upper = float(params.get('upper', 1.0))
                        lower = float(params.get('lower', 0.0))
                        sigma = float(params.get('sigma', 0.1))

                        import math
                        center = (upper + lower) / 2
                        # Gaussian centered at midpoint
                        result = math.exp(-((value - center) ** 2) / (2 * sigma ** 2))

                    else:
                        result = value

                    results.append(result)

                except (ValueError, TypeError) as e:
                    results.append(None)

            return Response({"results": results}, status=status.HTTP_200_OK)

        except Exception as e:
            log.exception("test modifier failed")
            return Response(
                {"error": str(e)},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR
            )



class PropertyScorerViewSet(viewsets.ModelViewSet):
    serializer_class = serializers.PropertyScorerSerializer

    def get_queryset(self):
        qs = models.PropertyScorer.objects.all().order_by("id")
        project_id = self.request.query_params.get("project_id")
        if project_id:
            qs = qs.filter(project_id=project_id)
        return qs

    @action(detail=False, methods=["get"], url_path="available-transforms")
    def available_transforms(self, request):
        """
        GET /reinvent/property-scorers/available-transforms/
        Returns all available REINVENT4 transform types with their metadata
        (description, parameters, defaults) and the per-property auto-defaults.
        """
        return Response({
            "transforms": models.REINVENT4_TRANSFORM_META,
            "property_defaults": models.REINVENT4_TRANSFORM_DEFAULTS,
        })

    @action(detail=False, methods=["post"], url_path="preview-transform")
    def preview_transform(self, request):
        """
        POST /reinvent/property-scorers/preview-transform/
        Body: {"transform": {...}, "x_min": float, "x_max": float, "n_points": int}
        Returns {"x": [...], "y": [...]} for plotting.
        """
        import numpy as np

        transform = request.data.get("transform")
        if not transform or not isinstance(transform, dict):
            return Response({"error": "transform dict required"}, status=status.HTTP_400_BAD_REQUEST)

        x_min = float(request.data.get("x_min", 0.0))
        x_max = float(request.data.get("x_max", 1.0))
        n_points = min(int(request.data.get("n_points", 200)), 1000)

        if x_min >= x_max:
            return Response({"error": "x_min must be < x_max"}, status=status.HTTP_400_BAD_REQUEST)

        x_values = list(np.linspace(x_min, x_max, n_points))
        try:
            y_values = models.compute_transform_values(transform, x_values)
        except Exception as e:
            log.exception("preview_transform compute failed")
            return Response({"error": str(e)}, status=status.HTTP_400_BAD_REQUEST)

        return Response({"x": x_values, "y": y_values})


class GenUIModelScorerViewSet(viewsets.ModelViewSet):
    serializer_class = serializers.GenUIModelScorerSerializer

    def get_queryset(self):
        qs = models.GenUIModelScorer.objects.all().order_by("id")
        project_id = self.request.query_params.get("project_id")
        if project_id:
            qs = qs.filter(project_id=project_id)
        return qs


class UnwantedSmartsScorerViewSet(viewsets.ModelViewSet):
    serializer_class = serializers.UnwantedSmartsScorerSerializer

    def get_queryset(self):
        qs = models.UnwantedSmartsScorer.objects.all().order_by("id")
        project_id = self.request.query_params.get("project_id")
        if project_id:
            qs = qs.filter(project_id=project_id)
        return qs



class ReinventEnvironmentViewSet(viewsets.ModelViewSet):
    queryset = models.ReinventEnvironment.objects.all()
    serializer_class = serializers.ReinventEnvironmentSerializer

    def get_queryset(self):
        qs = super().get_queryset().order_by("-id")
        project_id = self.request.query_params.get("project_id")
        if project_id:
            qs = qs.filter(project_id=project_id)
        return qs

    def perform_destroy(self, instance):
        # Delete dependent Reinvent runs and agents first, then the environment.
        # This is needed because the DB-level FK may still be RESTRICT on some
        # deployments where the migration did not update the constraint in-place.
        models.Reinvent.objects.filter(environment=instance).delete()
        models.ReinventAgent.objects.filter(environment=instance).delete()
        instance.delete()


class ReinventAgentTrainingViewSet(viewsets.ModelViewSet):
    queryset = models.ReinventAgentTraining.objects.all()
    serializer_class = serializers.ReinventAgentTrainingSerializer

    @action(detail=False, methods=["get"], url_path="learning-strategies")
    def learning_strategies(self, request):
        try:
            strategies = models.ReinventEnvironmentHelper.get_learning_strategies()
            return Response({"strategies": strategies})
        except Exception as e:
            log.exception("learning_strategies failed")
            return Response({"error": repr(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)


class ReinventAgentValidationViewSet(viewsets.ModelViewSet):
    queryset = models.ReinventAgentValidation.objects.all()
    serializer_class = serializers.ReinventAgentValidationSerializer


class ReinventAgentViewSet(viewsets.ModelViewSet):
    queryset = models.ReinventAgent.objects.all()
    serializer_class = serializers.ReinventAgentSerializer

    def get_queryset(self):
        qs = super().get_queryset()
        project_id = self.request.query_params.get("project_id")
        if project_id:
            qs = qs.filter(project_id=project_id)
        return qs


class ReinventStageViewSet(viewsets.ModelViewSet):
    queryset = models.ReinventStage.objects.all()
    serializer_class = serializers.ReinventStageSerializer

    def get_queryset(self):
        qs = super().get_queryset()
        generator_id = self.request.query_params.get("generator")
        if generator_id is not None:
            qs = qs.filter(generator_id=generator_id)
        return qs


class ReinventViewSet(viewsets.ModelViewSet):
    queryset = models.Reinvent.objects.order_by("-id")
    serializer_class = serializers.ReinventSerializer
    init_serializer_class = serializers.ReinventInitSerializer
    owner_relation = "project__owner"

    def get_queryset(self):
        qs = super().get_queryset()
        project_id = self.request.query_params.get("project_id")
        if project_id:
            qs = qs.filter(project_id=project_id)
        return qs

    def get_serializer_class(self):
        if self.action in {"create", "update", "partial_update"}:
            return self.init_serializer_class
        return self.serializer_class

    def create(self, request, *args, **kwargs):
        serializer = self.get_serializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        instance = serializer.save()

        build = bool(request.data.get("build", False))
        if build:
            device = _auto_detect_device(request.data.get("device", "auto"))
            eager = hasattr(settings, "CELERY_TASK_ALWAYS_EAGER") and settings.CELERY_TASK_ALWAYS_EAGER
            _, task_id = runTask(
                runReinventStagedLearning,
                instance=instance,
                eager=eager,
                args=(instance.id,),
                kwargs={"device": device},
            )
            out = serializers.ReinventSerializer(instance, context=self.get_serializer_context()).data
            out.update({"task_id": task_id, "device": device})
            return Response(out, status=status.HTTP_201_CREATED)

        out = serializers.ReinventSerializer(instance, context=self.get_serializer_context()).data
        return Response(out, status=status.HTTP_201_CREATED)

    @action(detail=True, methods=["post"], url_path="build-toml")
    def build_toml(self, request, pk=None):
        reinvent = self.get_object()
        device = _auto_detect_device(request.data.get("device", "auto"))
        try:
            toml_path = reinvent.build_staged_toml(device=device)
            return Response({"toml_path": toml_path}, status=status.HTTP_201_CREATED)
        except Exception as e:
            log.exception("build_toml failed for Reinvent pk=%s", pk)
            return Response({"error": repr(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

    @action(detail=True, methods=["post"], url_path="run-staged-learning")
    def run_staged_learning(self, request, pk=None):
        reinvent = self.get_object()
        device = _auto_detect_device(request.data.get("device", "auto"))
        try:
            eager = hasattr(settings, "CELERY_TASK_ALWAYS_EAGER") and settings.CELERY_TASK_ALWAYS_EAGER
            _, task_id = runTask(
                runReinventStagedLearning,
                instance=reinvent,
                eager=eager,
                args=(reinvent.id,),
                kwargs={"device": device},
            )
            return Response(
                {"task_id": task_id, "reinvent_id": reinvent.id, "device": device},
                status=status.HTTP_202_ACCEPTED,
            )
        except Exception as e:
            log.exception("run_staged_learning enqueue failed for Reinvent pk=%s", pk)
            return Response({"error": repr(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

    @action(detail=True, methods=["get"], url_path="download-results")
    def download_results(self, request, pk=None):
        """Download the CSV results from staged learning."""
        reinvent = self.get_object()
        try:
            mf = reinvent.get_results_file()
            needs_backfill = True
            if mf.file:
                try:
                    needs_backfill = mf.file.size == 0
                except Exception:
                    needs_backfill = True
            if needs_backfill:
                try:
                    csv_path = reinvent.agent.get_csv_path(generator=reinvent)
                    if os.path.isfile(csv_path):
                        with open(csv_path, "rb") as fh:
                            data = fh.read()
                        if data:
                            models._overwrite_filefield(mf, data, filename=os.path.basename(mf.file.name))
                except Exception:
                    pass

            if not mf.file:
                return Response({"error": "Results not available yet."}, status=status.HTTP_404_NOT_FOUND)
            try:
                size = mf.file.size
            except Exception:
                size = None
            if not size:
                return Response({"error": "Results not available yet."}, status=status.HTTP_404_NOT_FOUND)

            with mf.file.open("rb") as fh:
                csv_content = fh.read().decode("utf-8", errors="replace")

            return Response(csv_content, content_type="text/csv", status=status.HTTP_200_OK)
        except Exception as e:
            log.exception("download_results failed for Reinvent pk=%s", pk)
            return Response({"error": repr(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

    @action(detail=True, methods=["get"], url_path="csv-scores")
    def csv_scores(self, request, pk=None):
        """
        GET /reinvent/runs/{id}/csv-scores/
        Read ALL stage CSV files and return accurate score statistics + histogram bins.
        """
        import csv as csv_mod
        import glob
        import math

        reinvent = self.get_object()
        try:
            agent = reinvent.agent
            sl_dir = agent._sl_dir()
            prefix = agent._results_prefix(reinvent)
            legacy_prefix = getattr(agent.training, "summary_csv_prefix", None) or "reinvent"

            # Collect all numbered CSV files for this run
            def _find_csvs(pfx):
                pattern = os.path.join(sl_dir, f"{pfx}_*.csv")
                matches = glob.glob(pattern)
                plain = os.path.join(sl_dir, f"{pfx}.csv")
                if os.path.isfile(plain) and plain not in matches:
                    matches.append(plain)
                return matches

            csv_files = _find_csvs(prefix)
            if not csv_files:
                csv_files = _find_csvs(legacy_prefix)

            if not csv_files:
                return Response(
                    {"exists": False, "message": "No CSV result files found."},
                    status=status.HTTP_200_OK,
                )

            # Sort by stage number
            def _stage_num(p):
                stem = os.path.basename(p).rsplit(".", 1)[0]
                parts = stem.rsplit("_", 1)
                try:
                    return int(parts[-1])
                except (ValueError, IndexError):
                    return 0
            csv_files.sort(key=_stage_num)

            all_scores = []
            per_stage = []

            for csv_path in csv_files:
                if not os.path.isfile(csv_path):
                    continue
                stage_scores = []
                try:
                    with open(csv_path, "r", encoding="utf-8", newline="") as fh:
                        reader = csv_mod.reader(fh)
                        header_row = None
                        score_col = None
                        for row_num, row in enumerate(reader):
                            if len(row) < 2:
                                continue
                            if row_num == 0:
                                lower = [h.lower().strip() for h in row]
                                if any(h in ['smiles', 'step', 'agent', 'score'] for h in lower):
                                    header_row = lower
                                    for idx, h in enumerate(lower):
                                        if h in ('score', 'total_score', 'avg_score'):
                                            score_col = idx
                                            break
                                    continue
                            try:
                                if score_col is not None and score_col < len(row):
                                    val = float(row[score_col])
                                else:
                                    # Fallback: try second-to-last column
                                    val = float(row[-2])
                                if not math.isnan(val):
                                    stage_scores.append(val)
                            except (ValueError, IndexError):
                                continue
                except Exception as e:
                    log.warning(f"Could not parse CSV {csv_path}: {e}")
                    continue

                if stage_scores:
                    all_scores.extend(stage_scores)
                    stage_num = _stage_num(csv_path)
                    per_stage.append({
                        "stage": stage_num,
                        "file": os.path.basename(csv_path),
                        "count": len(stage_scores),
                        "mean": sum(stage_scores) / len(stage_scores),
                        "best": max(stage_scores),
                        "worst": min(stage_scores),
                        "median": sorted(stage_scores)[len(stage_scores) // 2],
                        "_scores": stage_scores,  # kept in memory for histogram, not sent to client
                    })

            if not all_scores:
                return Response(
                    {"exists": False, "message": "CSV files found but no scores parsed."},
                    status=status.HTTP_200_OK,
                )

            # Use the LAST stage for histogram and "latest" stats
            last_stage = per_stage[-1] if per_stage else None
            last_stage_scores = last_stage.get("_scores", all_scores) if last_stage else all_scores

            # Build histogram (20 bins) from LAST stage only
            hist_scores = last_stage_scores
            n_bins = 20
            min_s = min(hist_scores)
            max_s = max(hist_scores)
            bin_width = (max_s - min_s) / n_bins if max_s > min_s else 1.0
            if bin_width == 0:
                bin_width = 1.0
            bins = [0] * n_bins
            bin_edges = [min_s + i * bin_width for i in range(n_bins + 1)]
            for s in hist_scores:
                idx = int((s - min_s) / bin_width)
                if idx >= n_bins:
                    idx = n_bins - 1
                bins[idx] += 1


            payload = {
                "exists": True,
                "total_molecules": len(all_scores),
                "overall_mean": sum(all_scores) / len(all_scores),
                "overall_best": max(all_scores),
                "overall_worst": min(all_scores),
                "overall_median": sorted(all_scores)[len(all_scores) // 2],
                "latest_stage_mean": last_stage["mean"] if last_stage else None,
                "latest_stage_best": last_stage["best"] if last_stage else None,
                "stages": [{k: v for k, v in st.items() if k != "_scores"} for st in per_stage],
                "histogram": {
                    "bins": bins,
                    "edges": [round(e, 4) for e in bin_edges],
                    "n_bins": n_bins,
                    "stage": last_stage["stage"] if last_stage else None,
                },
            }
            return Response(payload, status=status.HTTP_200_OK)
        except Exception as e:
            log.exception("csv_scores failed for Reinvent pk=%s", pk)
            return Response({"error": repr(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

    @action(detail=True, methods=["get"], url_path="training-log")
    def training_log(self, request, pk=None):
        """
        GET /reinvent/runs/{id}/training-log/
        Returns the staged learning log content for this run.
        """
        try:
            run = self.get_object()
            agent = run.agent
            log_path = agent.get_rl_log_path(generator=run)

            if not os.path.isfile(log_path):
                return Response(
                    {"content": "", "exists": False, "message": "Staged learning log not available yet"},
                    status=status.HTTP_200_OK,
                )

            try:
                with open(log_path, "r", encoding="utf-8") as f:
                    content = f.read()
                return Response(
                    {"content": content, "exists": True, "path": log_path},
                    status=status.HTTP_200_OK,
                )
            except Exception as read_err:
                log.warning(f"Could not read staged learning log: {read_err}")
                return Response(
                    {"content": "", "exists": False, "message": f"Error reading log: {str(read_err)}"},
                    status=status.HTTP_200_OK,
                )
        except Exception as e:
            log.exception("training_log failed for Reinvent pk=%s", pk)
            return Response({"error": repr(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

    @action(detail=True, methods=["get"], url_path="training-tb")
    def training_tb(self, request, pk=None):
        """
        GET /reinvent/runs/{id}/training-tb/
        Returns summary statistics from TensorBoard logs for this staged learning run.
        REINVENT adds _0, _1, etc. suffixes to tb_logdir, so we search for matching folders.
        """
        try:
            run = self.get_object()
            agent = run.agent
            run_tag = agent._run_tag(run)
            sl_dir = agent._sl_dir()

            # REINVENT appends _0, _1, etc., so we need to find the folder
            tb_base = f"tb_logs_{run_tag}"
            tb_dir = None

            # Look for exact match first
            exact_path = os.path.join(sl_dir, tb_base)
            if os.path.isdir(exact_path):
                tb_dir = exact_path
            else:
                # Look for numbered variants: tb_logs_project5_run24_0, tb_logs_project5_run24_1, etc.
                import glob
                pattern = os.path.join(sl_dir, f"{tb_base}_*")
                matches = glob.glob(pattern)
                if matches:
                    # Use the most recent one
                    tb_dir = max(matches, key=os.path.getmtime)

            if not tb_dir or not os.path.isdir(tb_dir):
                return Response(
                    {"exists": False, "message": "TensorBoard log dir not found", "searched": tb_base},
                    status=status.HTTP_200_OK,
                )

            try:
                from tensorboard.backend.event_processing.event_accumulator import EventAccumulator
            except Exception as import_err:
                return Response(
                    {"exists": False, "message": f"TensorBoard reader unavailable: {str(import_err)}"},
                    status=status.HTTP_200_OK,
                )

            def _read_scalars(path: str):
                try:
                    ea = EventAccumulator(path)
                    ea.Reload()
                    tags = ea.Tags().get("scalars", [])
                    return ea, tags
                except Exception:
                    return None, []

            def _get_series(ea, tags, tag_candidates):
                if ea is None:
                    return None, []
                for tag in tag_candidates:
                    if tag in tags:
                        vals = ea.Scalars(tag)
                        if vals:
                            return tag, vals
                return None, []

            def _latest(vals):
                if not vals:
                    return None
                last = vals[-1]
                return {"step": last.step, "value": last.value, "wall_time": last.wall_time}

            def _best(vals):
                if not vals:
                    return None
                best = min(vals, key=lambda x: x.value)
                return {"step": best.step, "value": best.value, "wall_time": best.wall_time}

            def _series(vals, limit=500):
                if not vals:
                    return []
                tail = vals[-limit:] if limit and len(vals) > limit else vals
                return [{"step": v.step, "value": v.value, "wall_time": v.wall_time} for v in tail]

            # Root TB scalars
            ea_root, tags_root = _read_scalars(tb_dir)

            # Look for common staged learning metrics
            agent_nll_tag, agent_nll_vals = _get_series(ea_root, tags_root,
                ["Loss (likelihood averages)_agent NLL", "agent_nll", "Agent NLL"])
            prior_nll_tag, prior_nll_vals = _get_series(ea_root, tags_root,
                ["Loss (likelihood averages)_prior NLL", "prior_nll", "Prior NLL"])
            augmented_nll_tag, augmented_nll_vals = _get_series(ea_root, tags_root,
                ["Loss (likelihood averages)_augmented NLL", "augmented_nll", "Augmented NLL"])

            # Look for score metrics
            score_tag, score_vals = _get_series(ea_root, tags_root,
                ["score", "Score", "mean_score"])

            payload = {
                "exists": True,
                "tb_dir": tb_dir,
                "tags": tags_root,
                "agent_nll": {
                    "tag": agent_nll_tag,
                    "latest": _latest(agent_nll_vals),
                    "best": _best(agent_nll_vals),
                    "series": _series(agent_nll_vals),
                },
                "prior_nll": {
                    "tag": prior_nll_tag,
                    "latest": _latest(prior_nll_vals),
                    "best": _best(prior_nll_vals),
                    "series": _series(prior_nll_vals),
                },
                "augmented_nll": {
                    "tag": augmented_nll_tag,
                    "latest": _latest(augmented_nll_vals),
                    "best": _best(augmented_nll_vals),
                    "series": _series(augmented_nll_vals),
                },
                "score": {
                    "tag": score_tag,
                    "latest": _latest(score_vals),
                    "best": _best(score_vals) if score_vals else None,
                    "series": _series(score_vals),
                },
            }

            return Response(payload, status=status.HTTP_200_OK)
        except Exception as e:
            log.exception("training_tb failed for Reinvent pk=%s", pk)
            return Response({"error": repr(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

    @action(detail=True, methods=["get"], url_path="generated-molecules")
    def generated_molecules(self, request, pk=None):
        """
        GET /reinvent/runs/{id}/generated-molecules/
        Returns ALL generated molecule sets linked to this run,
        each with an explicit integer id and molecule count.
        """
        try:
            from genui.compounds.extensions.generated.models import GeneratedMolSet

            run = self.get_object()
            molsets = GeneratedMolSet.objects.filter(source=run).order_by("-id")

            if not molsets.exists():
                return Response(
                    {"exists": False, "message": "No generated molecules available yet", "molsets": []},
                    status=status.HTTP_200_OK,
                )

            result = []
            for ms in molsets:
                result.append({
                    "id": ms.pk,                        # always integer
                    "name": ms.name,
                    "description": ms.description or "",
                    "created": ms.created.isoformat() if ms.created else None,
                    "molecule_count": ms.molecules.count(),
                    "project": ms.project_id,
                })

            return Response(
                {"exists": True, "molsets": result},
                status=status.HTTP_200_OK,
            )
        except Exception as e:
            log.exception("generated_molecules failed for Reinvent pk=%s", pk)
            return Response({"error": repr(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

    @action(detail=True, methods=["delete"], url_path="generated-molecules/(?P<molset_pk>[0-9]+)")
    def delete_generated_molset(self, request, pk=None, molset_pk=None):
        """
        DELETE /reinvent/runs/{id}/generated-molecules/{molset_pk}/
        Deletes a generated molecule set belonging to this run.
        """
        try:
            from genui.compounds.extensions.generated.models import GeneratedMolSet

            run = self.get_object()
            try:
                ms = GeneratedMolSet.objects.get(pk=molset_pk, source=run)
            except GeneratedMolSet.DoesNotExist:
                return Response({"error": "Molecule set not found for this run."}, status=status.HTTP_404_NOT_FOUND)
            ms.delete()
            return Response(status=status.HTTP_204_NO_CONTENT)
        except Exception as e:
            log.exception("delete_generated_molset failed for Reinvent pk=%s molset_pk=%s", pk, molset_pk)
            return Response({"error": repr(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

    @action(detail=True, methods=["get"], url_path="generated-molecules/(?P<molset_pk>[0-9]+)/molecules")
    def molset_molecules(self, request, pk=None, molset_pk=None):
        """
        GET /reinvent/runs/{id}/generated-molecules/{molset_pk}/molecules/
        Returns paginated molecules from the molset with scores joined from the CSV.
        Query params: page (default 1), page_size (default 20), min_score (float)
        """
        import csv as csv_mod
        import glob as glob_mod
        import math as math_mod

        try:
            from genui.compounds.extensions.generated.models import GeneratedMolSet
            from genui.compounds.models import Molecule

            run = self.get_object()
            try:
                ms = GeneratedMolSet.objects.get(pk=molset_pk, source=run)
            except GeneratedMolSet.DoesNotExist:
                return Response({"error": "Molecule set not found for this run."}, status=status.HTTP_404_NOT_FOUND)

            # Build SMILES→score map from all CSV files for this run
            smiles_score: dict = {}
            try:
                agent = run.agent
                sl_dir = agent._sl_dir()
                prefix = agent._results_prefix(run)
                pattern = os.path.join(sl_dir, f"{prefix}_*.csv")
                csvs = glob_mod.glob(pattern)
                plain = os.path.join(sl_dir, f"{prefix}.csv")
                if os.path.isfile(plain):
                    csvs.append(plain)

                for csv_path in csvs:
                    if not os.path.isfile(csv_path):
                        continue
                    with open(csv_path, "r", encoding="utf-8", newline="") as fh:
                        reader = csv_mod.reader(fh)
                        score_col = smiles_col = None
                        for row_num, row in enumerate(reader):
                            if len(row) < 2:
                                continue
                            if row_num == 0:
                                lower = [h.lower().strip() for h in row]
                                for idx, h in enumerate(lower):
                                    if h == "smiles":
                                        smiles_col = idx
                                    if h in ("score", "total_score", "avg_score") and score_col is None:
                                        score_col = idx
                                if smiles_col is None:
                                    smiles_col = 4  # REINVENT default
                                if score_col is None:
                                    score_col = 3   # REINVENT default
                                continue
                            try:
                                smi = row[smiles_col].strip()
                                sc = float(row[score_col])
                                if smi and not math_mod.isnan(sc):
                                    # Keep best score per SMILES across all stages
                                    if smi not in smiles_score or sc > smiles_score[smi]:
                                        smiles_score[smi] = sc
                            except (IndexError, ValueError):
                                continue
            except Exception as csv_err:
                log.warning("molset_molecules: could not build score map: %s", csv_err)

            # Fetch all molecules from this molset
            mols_qs = Molecule.objects.filter(providers__id=ms.id).order_by("id")

            # Apply min_score filter
            min_score_param = request.query_params.get("min_score")
            min_score = None
            if min_score_param:
                try:
                    min_score = float(min_score_param)
                except ValueError:
                    pass

            if min_score is not None and smiles_score:
                # Filter in Python since scores aren't in DB
                all_mols = list(mols_qs)
                filtered = [m for m in all_mols if smiles_score.get(m.smiles, 0.0) >= min_score]
            else:
                filtered = list(mols_qs)

            total = len(filtered)

            # Paginate
            try:
                page = max(1, int(request.query_params.get("page", 1)))
                page_size = min(100, max(1, int(request.query_params.get("page_size", 20))))
            except ValueError:
                page, page_size = 1, 20

            start = (page - 1) * page_size
            page_mols = filtered[start:start + page_size]

            results = []
            for mol in page_mols:
                score = smiles_score.get(mol.smiles)
                results.append({
                    "id": mol.id,
                    "smiles": mol.smiles,
                    "inchi": mol.inchi,
                    "inchiKey": mol.inchiKey,
                    "score": round(score, 6) if score is not None else None,
                    "properties": mol.properties if hasattr(mol, 'properties') else {},
                })

            return Response({
                "count": total,
                "page": page,
                "page_size": page_size,
                "total_pages": math_mod.ceil(total / page_size) if page_size else 1,
                "results": results,
            }, status=status.HTTP_200_OK)

        except Exception as e:
            log.exception("molset_molecules failed for Reinvent pk=%s molset_pk=%s", pk, molset_pk)
            return Response({"error": repr(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)


    @action(detail=True, methods=["get"], url_path="tasks/all")
    def tasks_all(self, request, pk=None):
        """
        GET /reinvent/runs/{id}/tasks/all/
        Returns all tasks associated with this staged learning run.
        Uses getTasksAsDict() (the correct TaskShortcutsMixIn method) and
        cross-checks each task's live state via Celery AsyncResult so that
        actively running tasks are correctly reported even before
        django-celery-results has written the terminal state back.
        """
        from celery.result import AsyncResult

        try:
            run = self.get_object()
            tasks_dict = run.getTasksAsDict()  # {task_name: [{task_id, status, result, traceback}, ...]}

            serialized = []
            for task_name, task_list in tasks_dict.items():
                for t in task_list:
                    task_id = t.get("task_id")
                    db_status = str(t.get("status") or "").upper()

                    # Cross-check live Celery state — the DB record may lag
                    # (e.g. task is STARTED in the worker but DB still shows PENDING)
                    live_status = db_status
                    if task_id:
                        try:
                            ar = AsyncResult(task_id)
                            live = str(ar.state or "").upper()
                            # Prefer the live state when it carries more information
                            # Priority order: FAILURE > SUCCESS > STARTED/PROGRESS > PENDING > db value
                            priority = {"FAILURE": 5, "FAILED": 5, "REVOKED": 5,
                                        "SUCCESS": 4, "STARTED": 3, "PROGRESS": 3,
                                        "RECEIVED": 2, "PENDING": 1}
                            if priority.get(live, 0) >= priority.get(db_status, 0):
                                live_status = live
                        except Exception:
                            pass

                    serialized.append({
                        "id": task_id,
                        "status": live_status,
                        "name": task_name or "Staged Learning Task",
                        "result": t.get("result"),
                        "traceback": t.get("traceback"),
                    })

            return Response({"results": serialized}, status=status.HTTP_200_OK)
        except Exception as e:
            log.exception("tasks_all failed for Reinvent pk=%s", pk)
            return Response({"results": []}, status=status.HTTP_200_OK)


class ModelPerformanceReinventViewSet(viewsets.ReadOnlyModelViewSet):
    serializer_class = serializers.ModelPerformanceReinventSerializer

    def get_queryset(self):
        qs = models.ModelPerformanceReinvent.objects.all().order_by("created")

        agent_id = self.request.query_params.get("agent")
        if agent_id:
            qs = qs.filter(agent_id=agent_id)

        stage_index = self.request.query_params.get("stage_index")
        if stage_index is not None:
            try:
                qs = qs.filter(stage_index=int(stage_index))
            except (TypeError, ValueError):
                pass

        return qs


class ReinventModelFileViewSet(viewsets.ReadOnlyModelViewSet):
    serializer_class = ModelFileSerializer

    def get_queryset(self):
        qs = ModelFile.objects.all().order_by("-id")
        if getattr(self.request, "user", None) and not self.request.user.is_anonymous:
            qs = qs.filter(modelInstance__project__owner=self.request.user)

        project_id = self.request.query_params.get("project_id")
        if project_id:
            qs = qs.filter(modelInstance__project_id=project_id)

        return qs


class ReinventMetricsViewSet(MetricsViewSet):
    """
    Return only metrics valid for ReinventAgent mode validation.
    """
    def get_queryset(self):
        ret = super().get_queryset()
        modes = AlgorithmMode.objects.filter(name__in=("ReinventAgent",))
        return ret.filter(validModes__in=modes).distinct()
