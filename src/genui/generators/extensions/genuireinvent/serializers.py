# genui/generators/extensions/genuireinvent/serializers.py

from __future__ import annotations

from django.db.models import Q
from django.utils.translation import gettext_lazy as _
from rest_framework import serializers

from genui.compounds.models import MolSet
from genui.compounds.serializers import MolSetSerializer
from genui.models.models import ModelPerformanceMetric, Algorithm, AlgorithmMode, Model
from genui.models.serializers import (
    ModelSerializer,
    TrainingStrategyInitSerializer,
    TrainingStrategySerializer,
    ValidationStrategyInitSerializer,
    ValidationStrategySerializer,
)
from genui.projects.models import Project

from . import models


# =====================================================================
# 1) REINVENT NET (TL) STRATEGIES
# =====================================================================

class ReinventValidationStrategySerializer(ValidationStrategySerializer):
    class Meta:
        model = models.ReinventNetValidation
        fields = (
            ValidationStrategySerializer.Meta.fields
            + ("validSetSize", "split_method", "valid_fraction", "random_seed", "temporal_cutoff")
        )


class ReinventValidationStrategyInitSerializer(ValidationStrategyInitSerializer):
    # Optional explicit metrics
    metrics = serializers.PrimaryKeyRelatedField(
        many=True, queryset=ModelPerformanceMetric.objects.all(), required=False
    )

    class Meta:
        model = models.ReinventNetValidation
        fields = ReinventValidationStrategySerializer.Meta.fields


class ReinventTrainingStrategySerializer(TrainingStrategySerializer):
    class Meta(TrainingStrategySerializer.Meta):
        model = models.ReinventNetTraining

        _base_fields = tuple(getattr(TrainingStrategySerializer.Meta, "fields", ()))
        _extra_fields = (
            "epochs",
            "batch_size",
            "sample_batch_size",
            "save_every_n_epochs",
            "best_epoch",
            "best_valid_loss",
        )
        fields = _base_fields + _extra_fields

        _base_ro = tuple(getattr(TrainingStrategySerializer.Meta, "read_only_fields", ()))
        read_only_fields = _base_ro + (
            "best_epoch",
            "best_valid_loss",
            "started_at",
            "finished_at",
            "device",
        )


class ReinventTrainingStrategyInitSerializer(TrainingStrategyInitSerializer):
    class Meta(TrainingStrategyInitSerializer.Meta):
        model = models.ReinventNetTraining

        _base_fields = tuple(getattr(TrainingStrategyInitSerializer.Meta, "fields", ()))
        _extra_fields = (
            "epochs",
            "batch_size",
            "sample_batch_size",
            "save_every_n_epochs",
        )
        fields = _base_fields + _extra_fields


# =====================================================================
# 2) REINVENT NET SERIALIZERS
# =====================================================================

class ReinventNetSerializer(ModelSerializer):
    molset = MolSetSerializer(many=False, required=False, allow_null=True)

    # Override to avoid Model.trainingStrategy property throwing on duplicates.
    trainingStrategy = serializers.SerializerMethodField(read_only=True)
    validationStrategy = serializers.SerializerMethodField(read_only=True)

    parent = serializers.SerializerMethodField()
    bestEpoch = serializers.SerializerMethodField(read_only=True)
    bestValidLoss = serializers.SerializerMethodField(read_only=True)

    class Meta:
        model = models.ReinventNet
        # keep base fields, but don't expose performance here (usually read-only/handled elsewhere)
        fields = [f for f in ModelSerializer.Meta.fields if f != "performance"] + [
            "molset",
            "parent",
            "bestEpoch",
            "bestValidLoss",
        ]
        read_only_fields = [f for f in ModelSerializer.Meta.read_only_fields if f != "performance"] + [
            "bestEpoch",
            "bestValidLoss",
        ]

    def get_parent(self, obj):
        if obj.parent_id:
            return {"id": obj.parent_id, "name": getattr(obj.parent, "name", None)}
        return None

    def get_bestEpoch(self, obj):
        ts = obj.trainingStrategies.order_by("-id").first()
        return getattr(ts, "best_epoch", None) if ts else None

    def get_bestValidLoss(self, obj):
        ts = obj.trainingStrategies.order_by("-id").first()
        return getattr(ts, "best_valid_loss", None) if ts else None

    def get_trainingStrategy(self, obj):
        ts = obj.trainingStrategies.order_by("-id").first()
        return ReinventTrainingStrategySerializer(ts).data if ts else None

    def get_validationStrategy(self, obj):
        vs = obj.validationStrategies.order_by("-id").first()
        return ReinventValidationStrategySerializer(vs).data if vs else None


class ReinventNetInitSerializer(ReinventNetSerializer):
    molset = serializers.PrimaryKeyRelatedField(
        many=False, queryset=MolSet.objects.all(), required=False, allow_null=True
    )
    trainingStrategy = ReinventTrainingStrategyInitSerializer(many=False)
    validationStrategy = ReinventValidationStrategyInitSerializer(many=False, required=False)

    parent = serializers.PrimaryKeyRelatedField(
        many=False, queryset=models.ReinventNet.objects.all(), required=False, allow_null=True
    )

    class Meta:
        model = models.ReinventNet
        fields = ReinventNetSerializer.Meta.fields
        read_only_fields = ReinventNetSerializer.Meta.read_only_fields

    def create(self, validated_data, **kwargs):
        """
        Create ReinventNet + concrete training/validation strategies.

        IMPORTANT:
        - Do not remove trainingStrategy from validated_data before super().create(),
          because the base ModelSerializer may use it to resolve builder_id.
        """
        ts_data = dict(validated_data.get("trainingStrategy") or {})
        vs_data = dict(validated_data.get("validationStrategy") or {})
        parent = validated_data.get("parent")
        molset = validated_data.get("molset")

        instance = super().create(validated_data, molset=molset, **kwargs)

        if parent:
            instance.parent = parent
            instance.save(update_fields=["parent"])

        # Create concrete ReinventNetTraining
        if ts_data:
            trainingStrategy = models.ReinventNetTraining.objects.create(
                modelInstance=instance,
                algorithm=ts_data["algorithm"],
                mode=ts_data["mode"],
                epochs=ts_data.get(
                    "epochs",
                    models.ReinventNetTraining._meta.get_field("epochs").default,
                ),
                batch_size=ts_data.get(
                    "batch_size",
                    models.ReinventNetTraining._meta.get_field("batch_size").default,
                ),
                sample_batch_size=ts_data.get(
                    "sample_batch_size",
                    models.ReinventNetTraining._meta.get_field("sample_batch_size").default,
                ),
                save_every_n_epochs=ts_data.get(
                    "save_every_n_epochs",
                    models.ReinventNetTraining._meta.get_field("save_every_n_epochs").default,
                ),
            )
            self.saveParameters(trainingStrategy, ts_data)

        # Create concrete ReinventNetValidation (optional)
        if vs_data:
            validationStrategy = models.ReinventNetValidation.objects.create(
                modelInstance=instance,
                validSetSize=vs_data.get(
                    "validSetSize",
                    models.ReinventNetValidation._meta.get_field("validSetSize").default,
                ),
                split_method=vs_data.get(
                    "split_method",
                    models.ReinventNetValidation._meta.get_field("split_method").default,
                ),
                valid_fraction=vs_data.get(
                    "valid_fraction",
                    models.ReinventNetValidation._meta.get_field("valid_fraction").default,
                ),
                random_seed=vs_data.get(
                    "random_seed",
                    models.ReinventNetValidation._meta.get_field("random_seed").default,
                ),
                temporal_cutoff=vs_data.get("temporal_cutoff", None),
            )

            metrics = vs_data.get("metrics")
            if metrics:
                validationStrategy.metrics.set(metrics)
            else:
                # default metrics by algorithm/mode
                validationStrategy.metrics.set(
                    ModelPerformanceMetric.objects
                    .filter(validModes=ts_data.get("mode"))
                    .filter(
                        Q(validAlgorithms=ts_data.get("algorithm"))
                        | Q(validAlgorithms__isnull=True)
                    )
                    .distinct()
                )

            validationStrategy.save()

        return instance


# =====================================================================
# 3) RL / ENVIRONMENT CONFIG SERIALIZERS
# =====================================================================

class ReinventDiversityFilterSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.ReinventDiversityFilter
        fields = "__all__"


class ScoreModifierSerializer(serializers.Serializer):
    """
    Polymorphic serializer for ScoreModifier subclasses:
      - ClippedScore
      - SmoothHump
    """

    TYPE_CLIPPED = "ClippedScore"
    TYPE_HUMP = "SmoothHump"
    TYPE_CHOICES = (TYPE_CLIPPED, TYPE_HUMP)

    id = serializers.IntegerField(read_only=True)
    type = serializers.ChoiceField(choices=[(t, t) for t in TYPE_CHOICES])

    # DataSet fields (inherited)
    project = serializers.PrimaryKeyRelatedField(queryset=Project.objects.all())
    name = serializers.CharField(max_length=256)
    description = serializers.CharField(max_length=10000, required=False, allow_blank=True, allow_null=True)

    created = serializers.DateTimeField(read_only=True)
    updated = serializers.DateTimeField(read_only=True)

    # ---- ClippedScore fields ----
    upper = serializers.FloatField(required=False, allow_null=True)
    lower = serializers.FloatField(required=False, allow_null=True)
    high = serializers.FloatField(required=False, allow_null=True)
    low = serializers.FloatField(required=False, allow_null=True)
    smooth = serializers.BooleanField(required=False)

    # ---- SmoothHump fields ----
    sigma = serializers.FloatField(required=False, allow_null=True)

    def _downcast(self, obj):
        for attr in ("clippedscore", "smoothhump"):
            try:
                child = getattr(obj, attr)
            except Exception:
                child = None
            if child is not None:
                return child
        return obj

    def _model_for_type(self, t: str):
        if t == self.TYPE_CLIPPED:
            return models.ClippedScore
        if t == self.TYPE_HUMP:
            return models.SmoothHump
        raise serializers.ValidationError({"type": _("Unknown modifier type.")})

    def _base_repr(self, obj) -> dict:
        return {
            "id": obj.pk,
            "type": obj.__class__.__name__,
            "project": getattr(obj, "project_id", None),
            "name": getattr(obj, "name", None),
            "description": getattr(obj, "description", None),
            "created": getattr(obj, "created", None),
            "updated": getattr(obj, "updated", None),
        }

    def to_representation(self, obj):
        obj = self._downcast(obj)
        data = self._base_repr(obj)

        if isinstance(obj, models.ClippedScore):
            data.update(
                upper=obj.upper,
                lower=obj.lower,
                high=obj.high,
                low=obj.low,
                smooth=obj.smooth,
                sigma=None,
            )
        elif isinstance(obj, models.SmoothHump):
            data.update(
                upper=obj.upper,
                lower=obj.lower,
                high=None,
                low=None,
                smooth=None,
                sigma=obj.sigma,
            )
        else:
            data.update(
                upper=getattr(obj, "upper", None),
                lower=getattr(obj, "lower", None),
                high=getattr(obj, "high", None),
                low=getattr(obj, "low", None),
                smooth=getattr(obj, "smooth", None),
                sigma=getattr(obj, "sigma", None),
            )

        return data

    def validate(self, attrs):
        t = attrs.get("type")

        if t == self.TYPE_CLIPPED:
            upper = attrs.get("upper")
            lower = attrs.get("lower")
            if upper is None:
                raise serializers.ValidationError({"upper": _("This field is required for ClippedScore.")})
            if lower is not None and lower >= upper:
                raise serializers.ValidationError({"lower": _("lower must be < upper.")})
            return attrs

        if t == self.TYPE_HUMP:
            # model defaults handle missing values
            return attrs

        raise serializers.ValidationError({"type": _("Invalid modifier type.")})

    def create(self, validated_data):
        t = validated_data.pop("type")
        ModelCls = self._model_for_type(t)

        allowed = {f.name for f in ModelCls._meta.fields}
        kwargs = {k: v for k, v in validated_data.items() if k in allowed}

        return ModelCls.objects.create(**kwargs)

    def update(self, instance, validated_data):
        instance = self._downcast(instance)

        t = validated_data.pop("type", None)
        if t and t != instance.__class__.__name__:
            raise serializers.ValidationError(
                {"type": _("Changing modifier type is not supported. Create a new modifier instead.")}
            )

        for k, v in validated_data.items():
            if hasattr(instance, k):
                setattr(instance, k, v)

        instance.save()
        return instance



class PropertyScorerSerializer(serializers.ModelSerializer):
    # Read-only: the transform that will actually be written to the TOML.
    # Reflects transform_params if set, otherwise the auto-default for the property.
    effective_transform = serializers.SerializerMethodField(read_only=True)

    class Meta:
        model = models.PropertyScorer
        fields = "__all__"
        read_only_fields = ("id", "effective_transform")

    def get_effective_transform(self, obj):
        try:
            return obj.build_transform()
        except Exception:
            return None


class GenUIModelScorerSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.GenUIModelScorer
        fields = "__all__"
        read_only_fields = ("id",)


class UnwantedSmartsScorerSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.UnwantedSmartsScorer
        fields = "__all__"
        read_only_fields = ("id",)



class ReinventEnvironmentSerializer(serializers.ModelSerializer):
    prior_path = serializers.SerializerMethodField(read_only=True)
    agent_path = serializers.SerializerMethodField(read_only=True)
    prior_net_name = serializers.SerializerMethodField(read_only=True)
    agent_net_name = serializers.SerializerMethodField(read_only=True)
    diversity_filter_detail = serializers.SerializerMethodField(read_only=True)

    class Meta:
        model = models.ReinventEnvironment
        fields = "__all__"
        read_only_fields = ('id', 'created', 'updated', 'prior_path', 'agent_path',
                            'prior_net_name', 'agent_net_name', 'diversity_filter_detail')

    def create(self, validated_data):
        from django.utils import timezone
        # Ensure created and updated are set before calling model's __init__
        instance = models.ReinventEnvironment(**validated_data)
        if not instance.created:
            instance.created = timezone.now()
        if not instance.updated:
            instance.updated = timezone.now()
        instance.save()
        return instance

    def get_prior_path(self, obj):
        try:
            return obj.get_prior_path()
        except Exception:
            return None

    def get_agent_path(self, obj):
        try:
            return obj.get_agent_path()
        except Exception:
            return None

    def get_prior_net_name(self, obj):
        try:
            return obj.prior_net.name if obj.prior_net else None
        except Exception:
            return None

    def get_agent_net_name(self, obj):
        try:
            return obj.agent_net.name if obj.agent_net else None
        except Exception:
            return None

    def get_diversity_filter_detail(self, obj):
        try:
            df = obj.diversity_filter
            if not df:
                return None
            return {
                'id': df.id,
                'type': df.type,
                'bucket_size': df.bucket_size,
                'minscore': df.minscore,
                'minsimilarity': df.minsimilarity,
                'penalty_multiplier': df.penalty_multiplier,
            }
        except Exception:
            return None


# =====================================================================
# 4) RL AGENT CONFIG SERIALIZERS
# =====================================================================

class ReinventAgentTrainingSerializer(serializers.ModelSerializer):
    # Make TrainingStrategy base fields explicit so validation matches the DB constraints
    algorithm = serializers.PrimaryKeyRelatedField(
        queryset=Algorithm.objects.all(),
        required=True,
        allow_null=False
    )
    mode = serializers.PrimaryKeyRelatedField(
        queryset=AlgorithmMode.objects.all(),
        required=True,
        allow_null=False
    )
    modelInstance = serializers.PrimaryKeyRelatedField(
        queryset=Model.objects.all(),
        required=True,
        allow_null=False
    )
    # Explicit CharField for learning_type - will validate against model's choices at runtime
    learning_type = serializers.CharField(
        max_length=32,
        required=False,
        allow_blank=False,
        default="dap"
    )

    class Meta:
        model = models.ReinventAgentTraining
        fields = (
            "id",
            "algorithm",
            "mode",
            "modelInstance",
            "batch_size",
            "unique_sequences",
            "randomize_smiles",
            "tb_isim",
            "use_checkpoint",
            "purge_memories",
            "summary_csv_prefix",
            "learning_type",
            "sigma",
            "rate",
        )
        read_only_fields = ("id",)

    def to_internal_value(self, data):
        """Ensure learning_type defaults to 'dap' if not provided"""
        if not data.get('learning_type'):
            data = dict(data)  # Make a copy to avoid modifying original
            data['learning_type'] = 'dap'
        return super().to_internal_value(data)


class ReinventAgentValidationSerializer(serializers.ModelSerializer):
    # Make metrics optional with explicit field definition
    metrics = serializers.PrimaryKeyRelatedField(
        queryset=ModelPerformanceMetric.objects.all(),
        many=True,
        required=False,
        allow_empty=True
    )

    class Meta:
        model = models.ReinventAgentValidation
        fields = (
            "id",
            "modelInstance",
            "validate_every",
            "validation_dataset",
            "metrics",
        )
        read_only_fields = ("id",)

    def to_internal_value(self, data):
        """Ensure metrics defaults to empty list if not provided"""
        if not data.get('metrics'):
            data = dict(data)  # Make a copy to avoid modifying original
            data['metrics'] = []
        return super().to_internal_value(data)


class ReinventAgentSerializer(serializers.ModelSerializer):
    name = serializers.CharField(max_length=256, required=False, allow_blank=True)
    description = serializers.CharField(max_length=10000, required=False, allow_blank=True, default="")
    environment_name = serializers.SerializerMethodField(read_only=True)
    training_info = serializers.SerializerMethodField(read_only=True)
    validation_info = serializers.SerializerMethodField(read_only=True)

    class Meta:
        model = models.ReinventAgent
        fields = (
            "id",
            "name",
            "description",
            "environment",
            "environment_name",
            "training",
            "training_info",
            "validation",
            "validation_info",
            "output_model",
            "tb_logdir",
            "json_out_config",
        )
        read_only_fields = ("id", "output_model", "environment_name", "training_info", "validation_info")

    def get_environment_name(self, obj):
        env = getattr(obj, "environment", None)
        return env.name if env else None

    def get_training_info(self, obj):
        t = getattr(obj, "training", None)
        if not t:
            return None
        return {
            "id": t.id,
            "batch_size": getattr(t, "batch_size", None),
            "learning_type": getattr(t, "learning_type", None),
            "sigma": getattr(t, "sigma", None),
            "rate": getattr(t, "rate", None),
        }

    def get_validation_info(self, obj):
        v = getattr(obj, "validation", None)
        if not v:
            return None
        return {
            "id": v.id,
            "validate_every": getattr(v, "validate_every", None),
        }

    def create(self, validated_data):
        """
        Create ReinventAgent with proper project assignment.
        The project is extracted from the related environment/training/validation objects.
        """
        from genui.models.models import ModelBuilder

        # Get the environment to extract project
        environment = validated_data.get('environment')

        if not environment:
            raise serializers.ValidationError({"environment": "Environment is required to determine the project."})

        # Extract project from environment
        project = environment.project

        # Add project to validated_data if not present
        validated_data['project'] = project

        # Generate a default name if not provided
        if 'name' not in validated_data or not validated_data.get('name'):
            validated_data['name'] = f"Agent_{environment.id}"

        # Set builder to ReinventAgent if not provided (required field)
        if 'builder' not in validated_data or not validated_data.get('builder'):
            try:
                # Try to get the ReinventAgent builder
                builder = ModelBuilder.objects.get(name='ReinventAgent')
                validated_data['builder'] = builder
            except ModelBuilder.DoesNotExist:
                # If not found, get the first available builder
                builder = ModelBuilder.objects.first()
                if builder:
                    validated_data['builder'] = builder
                else:
                    raise serializers.ValidationError(
                        {"builder": "No ModelBuilder found in database. ReinventAgent builder must be registered."}
                    )

        # Call parent create method
        return super().create(validated_data)


# =====================================================================
# 5) STAGED LEARNING / RL RUN SERIALIZERS
# =====================================================================

class ReinventStageSerializer(serializers.ModelSerializer):
    resolved_checkpoint_path = serializers.SerializerMethodField(read_only=True)

    class Meta:
        model = models.ReinventStage
        fields = (
            "id",
            "generator",
            "order",
            "termination_type",
            "max_score",
            "min_steps",
            "max_steps",
            "scoring_source",
            "aggregation_type",
            "scoring_file",
            "chkpt_net",
            "chkpt_file",
            "property_scorers",
            "model_scorers",
            "smarts_scorers",
            "resolved_checkpoint_path",
        )
        read_only_fields = ("id", "scoring_file", "chkpt_file", "resolved_checkpoint_path")

    def get_resolved_checkpoint_path(self, obj):
        try:
            return obj.resolve_checkpoint_path(default_path="")
        except Exception:
            return None


class ReinventSerializer(serializers.ModelSerializer):
    agent_name = serializers.SerializerMethodField(read_only=True)
    environment_name = serializers.SerializerMethodField(read_only=True)

    class Meta:
        model = models.Reinvent
        fields = (
            "id",
            "name",
            "description",
            "created",
            "agent",
            "agent_name",
            "environment",
            "environment_name",
        )
        read_only_fields = ("id", "agent_name", "environment_name", "created")

    def get_agent_name(self, obj):
        try:
            return obj.agent.name if obj.agent else None
        except Exception:
            return None

    def get_environment_name(self, obj):
        try:
            return obj.environment.name if obj.environment else None
        except Exception:
            return None


class ReinventInlineStageSerializer(serializers.Serializer):
    """Lightweight serializer for stages embedded in run creation payload."""
    order = serializers.IntegerField(default=1)
    termination_type = serializers.CharField(default="simple", max_length=64)
    max_score = serializers.FloatField(default=1.0)
    min_steps = serializers.IntegerField(default=1)
    max_steps = serializers.IntegerField(default=100)
    scoring_source = serializers.CharField(default="inline", max_length=16)
    aggregation_type = serializers.CharField(default="geometric_mean", max_length=64, required=False)
    property_scorers = serializers.PrimaryKeyRelatedField(
        many=True, queryset=models.PropertyScorer.objects.all(), required=False
    )
    model_scorers = serializers.PrimaryKeyRelatedField(
        many=True, queryset=models.GenUIModelScorer.objects.all(), required=False
    )
    smarts_scorers = serializers.PrimaryKeyRelatedField(
        many=True, queryset=models.UnwantedSmartsScorer.objects.all(), required=False
    )


class ReinventInitSerializer(ReinventSerializer):
    # write-through fields (stored on linked agent)
    tb_logdir = serializers.CharField(required=False, allow_blank=True, write_only=True)
    json_out_config = serializers.CharField(required=False, allow_blank=True, write_only=True)

    # mirror TL workflow: allow build=true on POST
    build = serializers.BooleanField(required=False, default=False, write_only=True)
    device = serializers.CharField(required=False, default="cuda:0", write_only=True)

    # Bad SMARTS penalty weight (0–1)
    bad_smarts_weight = serializers.FloatField(required=False, default=1.0, write_only=True)

    # Inline stages — created together with the run
    stages = ReinventInlineStageSerializer(many=True, required=False, write_only=True)

    # Make environment optional since it comes from the agent
    environment = serializers.PrimaryKeyRelatedField(
        queryset=models.ReinventEnvironment.objects.all(),
        required=False,
        allow_null=True
    )

    class Meta(ReinventSerializer.Meta):
        fields = (
            "id",
            "name",
            "description",
            "agent",
            "environment",
            "tb_logdir",
            "json_out_config",
            "build",
            "device",
            "bad_smarts_weight",
            "stages",
        )
        read_only_fields = ("id",)

    def create(self, validated_data):
        from django.db import transaction

        tb_logdir = validated_data.pop("tb_logdir", None)
        json_out_config = validated_data.pop("json_out_config", None)
        validated_data.pop("build", None)
        device = validated_data.pop("device", None)
        stages_data = validated_data.pop("stages", [])
        # bad_smarts_weight stays in validated_data — it's a model field

        # Extract project and environment from agent
        agent = validated_data.get('agent')
        if agent:
            # Get project from agent
            validated_data['project'] = agent.project
            # Auto-populate environment from agent if not provided
            if 'environment' not in validated_data or validated_data.get('environment') is None:
                validated_data['environment'] = agent.environment
        else:
            raise serializers.ValidationError({"agent": "Agent is required to determine the project."})

        with transaction.atomic():
            instance = super().create(validated_data)

            # write-through to agent
            agent_instance = getattr(instance, "agent", None)
            if agent_instance and (tb_logdir is not None or json_out_config is not None):
                fields = []
                if tb_logdir is not None:
                    agent_instance.tb_logdir = tb_logdir
                    fields.append("tb_logdir")
                if json_out_config is not None:
                    agent_instance.json_out_config = json_out_config
                    fields.append("json_out_config")
                if fields:
                    agent_instance.save(update_fields=fields)

            # Create inline stages with scorer assignments
            for stage_data in stages_data:
                prop_scorers = stage_data.pop("property_scorers", [])
                model_scorers = stage_data.pop("model_scorers", [])
                smarts_scorers = stage_data.pop("smarts_scorers", [])

                stage = models.ReinventStage.objects.create(
                    generator=instance,
                    **stage_data,
                )
                if prop_scorers:
                    stage.property_scorers.set(prop_scorers)
                if model_scorers:
                    stage.model_scorers.set(model_scorers)
                if smarts_scorers:
                    stage.smarts_scorers.set(smarts_scorers)

        return instance


# =====================================================================
# 6) PERFORMANCE LOGGING
# =====================================================================

class ModelPerformanceReinventSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.ModelPerformanceReinvent
        fields = "__all__"

