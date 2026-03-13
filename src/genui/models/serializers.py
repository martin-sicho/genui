import inspect
from abc import abstractmethod

from rest_framework import serializers
from rest_polymorphic.serializers import PolymorphicSerializer

from genui.utils.serializers import GenericModelSerializerMixIn
from genui.models.models import (ModelFileFormat, ModelBuilder, Model, PARAM_VALUE_CTYPE_TO_MODEL_MAP, ModelParameter, \
                                 Algorithm, TrainingStrategy, ModelFile, BasicValidationStrategy, ValidationStrategy, \
                                 AlgorithmMode, ModelParameterValue, ModelPerformance,
                                 HyperparameterOptimizationStrategy, GridSearchOptimization, OptunaOptimization)
from genui.models import models
from genui.projects.models import Project


class ModelFileFormatSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = ModelFileFormat
        fields = ('id', 'fileExtension', 'description')


class ModelParameterSerializer(serializers.HyperlinkedModelSerializer):
    defaultValue = serializers.SerializerMethodField(read_only=True)

    class Meta:
        model = ModelParameter
        fields = ('id', 'name', 'contentType', 'defaultValue')
        read_only_fields = ('defaultValue',)

    def get_defaultValue(self, obj):
        return obj.defaultValue.value


class AlgorithmModeSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = AlgorithmMode
        fields = ('id', 'name',)


class AlgorithmSerializer(serializers.HyperlinkedModelSerializer):
    fileFormats = ModelFileFormatSerializer(many=True)
    parameters = ModelParameterSerializer(many=True)
    validModes = AlgorithmModeSerializer(many=True)

    class Meta:
        model = Algorithm
        fields = ('id', 'name', 'fileFormats', 'parameters', 'validModes')


class ModelPerformanceSerializer(GenericModelSerializerMixIn, serializers.HyperlinkedModelSerializer):
    className = GenericModelSerializerMixIn.className
    extraArgs = GenericModelSerializerMixIn.extraArgs
    metric = serializers.CharField()

    class Meta:
        model = ModelPerformance
        fields = ('id', 'value', 'metric', 'className', 'extraArgs')


class ModelParameterValueSerializer(serializers.HyperlinkedModelSerializer):
    parameter = ModelParameterSerializer(many=False)
    value = serializers.CharField()

    class Meta:
        model = ModelParameterValue
        fields = ('id', 'parameter', 'value')


class HyperparameterOptimizationStrategySerializer(serializers.HyperlinkedModelSerializer):
    metric = serializers.CharField()
    scoreAggregation = serializers.CharField()
    searchSpace = serializers.JSONField()

    class Meta:
        model = HyperparameterOptimizationStrategy
        fields = ('searchSpace', 'scoreAggregation', 'trainingStrategy', 'metric')


class HyperparameterOptimizationStrategyInitSerializer(HyperparameterOptimizationStrategySerializer):
    metric = serializers.CharField()
    scoreAggregation = serializers.CharField()
    searchSpace = serializers.JSONField()

    allowed_types = []

    class Meta:
        model = HyperparameterOptimizationStrategy
        fields = tuple(x for x in HyperparameterOptimizationStrategySerializer.Meta.fields if x != 'trainingStrategy')

    def validate_searchSpace_param(self, value):
        if "name" not in value or "type" not in value or "value" not in value:
            raise serializers.ValidationError(f"Invalid parameter in search space: {value}."
                                              f" Must contain 'name', 'type' and 'value'.")
        if value["type"] not in self.allowed_types:
            raise serializers.ValidationError(f"Invalid type for parameter {value['name']}: {value['type']}. "
                                              f" Can be one of: {self.allowed_types}.")
        return value["name"], value["type"], value["value"]


class GridSearchOptimizationInitSerializer(HyperparameterOptimizationStrategyInitSerializer):
    class Meta:
        model = GridSearchOptimization
        fields = HyperparameterOptimizationStrategyInitSerializer.Meta.fields

    allowed_types = ["range", "sequence"]

    def validate_searchSpace(self, value):
        for param in value:
            param_name, param_type, param_value = self.validate_searchSpace_param(param)
            if param_type == "sequence":
                if not isinstance(param_value, list):
                    raise serializers.ValidationError(
                        f"Invalid value of parameter {param_name} in search space: {param_value}."
                        f" Must be a list of values")
                if len(param_value) != len(set(param_value)):
                    raise serializers.ValidationError(
                        f"Invalid value of parameter {param_name} in search space: {param_value}."
                        f" List of values must not contain duplicates")
            elif param_type == "range":
                if not 2 <= len(param_value) <= 3:
                    raise serializers.ValidationError(
                        f"Invalid value of parameter {param_name} in search space: {param_value}."
                        f" Must be a list of two ot three values")
                if len(param_value) == 2:
                    start, stop = param_value
                else:
                    start, stop, step = param_value
                    if not isinstance(step, int) or step <= 0:
                        raise serializers.ValidationError(
                            f"Invalid value of parameter {param_name} in search space: {param_value}."
                            f" Step must be positive non-zero integer")
                if not isinstance(start, int) or not isinstance(stop, int):
                    raise serializers.ValidationError(
                        f"Invalid value of parameter {param_name} in search space: {param_value}."
                        f" Both start and stop must be integers")
                if start >= stop:
                    raise serializers.ValidationError(
                        f"Invalid value of parameter {param_name} in search space: {param_value}."
                        f" Stop must be greater than start.")
        return value


class GridSearchOptimizationSerializer(GridSearchOptimizationInitSerializer):
    class Meta:
        model = GridSearchOptimization
        fields = GridSearchOptimizationInitSerializer.Meta.fields


class OptunaOptimizationInitSerializer(HyperparameterOptimizationStrategyInitSerializer):
    nTrials = serializers.IntegerField(min_value=1)

    allowed_types = ["int", "float", "categorical"]

    class Meta:
        model = OptunaOptimization
        fields = HyperparameterOptimizationStrategyInitSerializer.Meta.fields + ('nTrials',)

    def validate_searchSpace(self, value):
        for param in value:
            param_name, param_type, param_value = self.validate_searchSpace_param(param)
            if param_type in ["int", "float"]:
                min_, max_ = param_value
                type_ = {"int": int, "float": float}[param_type]
                min_ = type_(min_)
                max_ = type_(max_)
                if min_ >= max_:
                    raise serializers.ValidationError(f"{param_name} in search space:"
                                                      " First value must be less than the second.")
            elif param_type == "categorical":
                if len(param_value) > len(set(param_value)):
                    raise serializers.ValidationError(f"{param_name} in search space:"
                                                      f" List of values must not contain duplicates")
        return value


class OptunaOptimizationSerializer(OptunaOptimizationInitSerializer):
    nTrials = serializers.IntegerField(min_value=1)

    class Meta:
        model = OptunaOptimization
        fields = OptunaOptimizationInitSerializer.Meta.fields


class HyperparameterOptimizationStrategyPolymorphicSerializer(PolymorphicSerializer):
    model_serializer_mapping = {
        GridSearchOptimization: GridSearchOptimizationSerializer,
        OptunaOptimization: OptunaOptimizationSerializer,
    }


class HyperparameterOptimizationStrategyPolymorphicInitSerializer(PolymorphicSerializer):
    model_serializer_mapping = {
        HyperparameterOptimizationStrategy: HyperparameterOptimizationStrategyInitSerializer,
        GridSearchOptimization: GridSearchOptimizationInitSerializer,
        OptunaOptimization: OptunaOptimizationInitSerializer,
    }


class ValidationStrategySerializer(serializers.HyperlinkedModelSerializer):
    metrics = serializers.CharField()

    class Meta:
        model = ValidationStrategy
        fields = ("metrics", "trainingStrategy")


class ValidationStrategyInitSerializer(ValidationStrategySerializer):
    metrics = serializers.CharField()

    class Meta:
        model = ValidationStrategy
        fields = tuple(x for x in ValidationStrategySerializer.Meta.fields if x != 'trainingStrategy')


class BasicValidationStrategyInitSerializer(ValidationStrategyInitSerializer):
    metrics = serializers.ListSerializer(child=serializers.CharField())
    cvFolds = serializers.IntegerField(min_value=0)
    dataSplit = serializers.JSONField()

    # TODO: check if correct metrics are used with the correct algorithm

    class Meta:
        model = BasicValidationStrategy
        fields = ValidationStrategyInitSerializer.Meta.fields + ('cvFolds', 'dataSplit')

    def validate_dataSplit(self, value):
        def validate_split_config(config):
            if not isinstance(config, dict):
                raise serializers.ValidationError("Split configuration must be a dictionary")

            if "name" not in config:
                raise serializers.ValidationError("Missing required 'name' key")

            try:
                module_path, class_name = config["name"].rsplit('.', 1)
                module = __import__(module_path, fromlist=[class_name])
                split_class = getattr(module, class_name)
                if inspect.isabstract(split_class):
                    raise serializers.ValidationError(f"Split class '{class_name}' cannot be abstract")
            except (ImportError, AttributeError) as e:
                raise serializers.ValidationError(f"Invalid split function: {str(e)}")

            valid_params = {
                param for param in config.keys()
                if param != "name" and param != "method"
            }

            for param in valid_params:
                if isinstance(config[param], dict):
                    validate_split_config(config[param])

            sig = inspect.signature(split_class.__init__)
            allowed_params = {
                param.name for param in sig.parameters.values()
                if param.name != 'self'
            }

            invalid_params = valid_params - allowed_params
            if invalid_params:
                raise serializers.ValidationError(
                    f"Invalid parameters for {class_name}: {', '.join(invalid_params)}"
                )

            return True

        try:
            validate_split_config(value)
        except Exception as e:
            raise serializers.ValidationError(str(e))

        return value


class BasicValidationStrategySerializer(BasicValidationStrategyInitSerializer):
    metrics = serializers.ListSerializer(child=serializers.CharField())


class ValidationStrategyPolymorphicSerializer(PolymorphicSerializer):
    model_serializer_mapping = {
        BasicValidationStrategy: BasicValidationStrategySerializer,
        # Add other subclasses and their serializers here
    }


class ValidationStrategyPolymorphicInitSerializer(PolymorphicSerializer):
    model_serializer_mapping = {
        ValidationStrategy: ValidationStrategyInitSerializer,
        BasicValidationStrategy: BasicValidationStrategyInitSerializer,
    }


class TrainingStrategySerializer(serializers.HyperlinkedModelSerializer):
    algorithm = AlgorithmSerializer(many=False)
    parameters = ModelParameterValueSerializer(many=True)
    mode = AlgorithmModeSerializer(many=False)
    validationStrategies = ValidationStrategyPolymorphicSerializer(many=True, read_only=True)
    hyperParamOptStrategies = HyperparameterOptimizationStrategyPolymorphicSerializer(many=True, read_only=True)

    class Meta:
        model = TrainingStrategy
        fields = ('id', 'algorithm', 'mode', 'parameters', 'validationStrategies', 'hyperParamOptStrategies')


class TrainingStrategyInitSerializer(TrainingStrategySerializer):
    algorithm = serializers.PrimaryKeyRelatedField(many=False, queryset=Algorithm.objects.all())
    parameters = serializers.DictField(allow_empty=True, child=serializers.CharField(), required=False)
    mode = serializers.PrimaryKeyRelatedField(many=False, queryset=AlgorithmMode.objects.all())
    validationStrategies = ValidationStrategyPolymorphicInitSerializer(many=True, required=False)
    hyperParamOptStrategies = HyperparameterOptimizationStrategyPolymorphicInitSerializer(many=True, required=False)

    class Meta:
        model = TrainingStrategy
        fields = TrainingStrategySerializer.Meta.fields

    def create(self, validated_data, **kwargs):
        instance = super().create(validated_data, **kwargs)

        if 'validationStrategies' in validated_data:
            strat_data = validated_data['validationStrategies']
            for vs_data in strat_data:
                validationStrategy = BasicValidationStrategy.objects.create(
                    trainingStrategy=instance,
                    cvFolds=vs_data['cvFolds'],
                    dataSplit=vs_data['dataSplit'],
                )
                validationStrategy.metrics.set(vs_data['metrics'])
                validationStrategy.save()

        if 'hyperParamOptStrategies' in validated_data:
            hypo_data = validated_data['hyperParamOptStrategies'][0]
            hyperparam_opt_strategy_class = getattr(models, hypo_data["resourcetype"])
            hypo_data.pop("resourcetype")
            hyperparam_opt_strategy = hyperparam_opt_strategy_class.objects.create(
                trainingStrategy=instance,
                **hypo_data
            )
            hyperparam_opt_strategy.save()

        return instance


class ModelFileSerializer(serializers.HyperlinkedModelSerializer):
    model = serializers.PrimaryKeyRelatedField(many=False, required=False, queryset=Model.objects.all())
    format = ModelFileFormatSerializer(many=False, read_only=True)

    class Meta:
        model = ModelFile
        fields = ('id', 'file', 'format', 'kind', 'model', 'note')
        read_only_fields = ('id', 'format',)

    def create(self, validated_data):
        return ModelFile.create(
            validated_data['model'],
            validated_data['file'].name,
            validated_data['file'],
            validated_data['kind'] if 'kind' in validated_data else ModelFile.AUXILIARY,
            validated_data['note'] if 'note' in validated_data else ''
        )


class ModelSerializer(serializers.HyperlinkedModelSerializer):
    project = serializers.PrimaryKeyRelatedField(many=False, queryset=Project.objects.all())
    trainingStrategy = TrainingStrategySerializer(many=False)
    build = serializers.BooleanField(default=True)
    taskID = serializers.UUIDField(required=False, read_only=True, allow_null=True)
    modelFile = ModelFileSerializer(many=False, read_only=True, allow_null=True, required=False)

    def __init__(self, *args, builder_class=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.builder_class = self.instance.builder.name if self.instance and isinstance(self.instance,
                                                                                        Model) else builder_class

    def is_valid(self, *, raise_exception=False):
        ret = super().is_valid(raise_exception=raise_exception)
        if ("trainingStrategy" in self.validated_data
                and "hyperParamOptStrategies" in self.validated_data["trainingStrategy"]
                and len(self.validated_data["trainingStrategy"]["hyperParamOptStrategies"]) > 0):
            val_data = self.validated_data["trainingStrategy"]["validationStrategies"]
            if len(val_data) > 1:
                raise serializers.ValidationError(
                    "Only one validation strategy is allowed when using hyperparameter optimization.")
        return ret

    @staticmethod
    def saveParameters(strat_instance: TrainingStrategy, strat_data):
        if 'parameters' not in strat_data:
            return

        alg_name = strat_data['algorithm'].name
        alg = Algorithm.objects.get(name=alg_name)

        for param in alg.parameters.all():
            if param.name not in strat_data['parameters']:
                strat_data['parameters'][param.name] = str(param.defaultValue.value)

        for param_name in strat_data['parameters']:
            parameter = ModelParameter.objects.get(
                name=param_name
                , algorithm__name=alg_name
            )
            value_class = PARAM_VALUE_CTYPE_TO_MODEL_MAP[parameter.contentType]
            parameter_value = value_class(
                parameter=parameter
                , strategy=strat_instance
                , value=value_class.parseValue(strat_data['parameters'][param_name]))
            parameter_value.save()

    class Meta:
        model = Model
        fields = (
            'id', 'name', 'description', 'created', 'updated', 'project', 'trainingStrategy', 'modelFile', 'build',
            'taskID')
        read_only_fields = ('id', 'created', 'updated', 'modelFile', 'taskID')

    def useBuilder(self, builder_class):
        self.builder_class = builder_class

    def create(self, validated_data, **kwargs):
        instance = self.Meta.model.objects.create(
            name=validated_data['name'],
            description=validated_data['description'] if 'description' in validated_data else '',
            project=validated_data['project'],
            builder=ModelBuilder.objects.get_or_create(
                name=self.builder_class if type(self.builder_class) == str else self.builder_class.__name__
            )[0],
            **kwargs
        )
        instance.build = validated_data["build"]
        return instance
