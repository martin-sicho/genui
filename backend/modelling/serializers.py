"""
serializers

Created by: Martin Sicho
On: 24-01-20, 14:44
"""
from rest_framework import serializers
import uuid

import modelling.models
from commons.serializers import GenericModelSerializerMixIn
from projects.models import Project


class ModelPerformanceMetricSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = modelling.models.ModelPerformanceMetric
        fields = ('id', 'name', 'description')


class ModelPerformanceSerializer(GenericModelSerializerMixIn, serializers.HyperlinkedModelSerializer):
    className = GenericModelSerializerMixIn.className
    extraArgs = GenericModelSerializerMixIn.extraArgs
    metric = ModelPerformanceMetricSerializer(many=False)

    class Meta:
        model = modelling.models.ModelPerformance
        fields = ('id', 'value', 'metric', 'className', 'extraArgs')


class ModelFileFormatSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = modelling.models.ModelFileFormat
        fields = ('id', 'fileExtension', 'description')


class ModelParameterSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = modelling.models.ModelParameter
        fields = ('id', 'name', 'contentType')


class AlgorithmModeSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = modelling.models.AlgorithmMode
        fields = ('id', 'name',)


class AlgorithmSerializer(serializers.HyperlinkedModelSerializer):
    fileFormats = ModelFileFormatSerializer(many=True)
    parameters = ModelParameterSerializer(many=True)
    validModes = AlgorithmModeSerializer(many=True)

    class Meta:
        model = modelling.models.Algorithm
        fields = ('id', 'name', 'fileFormats', 'parameters', 'validModes')


class ModelParameterValueSerializer(serializers.HyperlinkedModelSerializer):
    parameter = ModelParameterSerializer(many=False)
    value = serializers.CharField()

    class Meta:
        model = modelling.models.ModelParameterValue
        fields = ('id','parameter', 'value')


class TrainingStrategySerializer(serializers.HyperlinkedModelSerializer):
    algorithm = AlgorithmSerializer(many=False)
    parameters = ModelParameterValueSerializer(many=True)
    mode = AlgorithmModeSerializer(many=False)

    class Meta:
        model = modelling.models.TrainingStrategy
        fields = ('algorithm', 'mode', 'parameters')

class TrainingStrategyInitSerializer(TrainingStrategySerializer):
    algorithm = serializers.PrimaryKeyRelatedField(many=False, queryset=modelling.models.Algorithm.objects.all())
    parameters = serializers.DictField(allow_empty=True, child=serializers.CharField(), required=False)
    mode = serializers.PrimaryKeyRelatedField(many=False, queryset=modelling.models.AlgorithmMode.objects.all())

    class Meta:
        model = modelling.models.TrainingStrategy
        fields = TrainingStrategySerializer.Meta.fields

class ValidationStrategySerializer(serializers.HyperlinkedModelSerializer):
    metrics = ModelPerformanceMetricSerializer(many=True)

    class Meta:
        model = modelling.models.ValidationStrategy
        fields = ("metrics",)

class ValidationStrategyInitSerializer(ValidationStrategySerializer):
    metrics = serializers.PrimaryKeyRelatedField(many=True, queryset=modelling.models.ModelPerformanceMetric.objects.all())

    class Meta:
        model = modelling.models.ValidationStrategy
        fields = ValidationStrategySerializer.Meta.fields

class BasicValidationStrategyInitSerializer(ValidationStrategySerializer):
    metrics = serializers.PrimaryKeyRelatedField(many=True, queryset=modelling.models.ModelPerformanceMetric.objects.all())
    cvFolds = serializers.IntegerField(min_value=0)
    validSetSize = serializers.FloatField(min_value=0)

    class Meta:
        model = modelling.models.BasicValidationStrategy
        fields = ValidationStrategySerializer.Meta.fields + ('cvFolds', 'validSetSize')

class ModelFileSerializer(serializers.HyperlinkedModelSerializer):
    model = serializers.PrimaryKeyRelatedField(required=False, queryset=modelling.models.Model.objects.all())
    format = ModelFileFormatSerializer(many=False, read_only=True)

    class InvalidFileFormatError(Exception):

        def __init__(self, msg):
            super().__init__(msg)

    class Meta:
        model = modelling.models.ModelFile
        fields = ('file', 'format', 'kind', 'model')
        read_only_fields = ('format',)

    def create(self, validated_data):
        ModelFile = modelling.models.ModelFile
        model = validated_data['model']
        if validated_data['kind'] == ModelFile.MAIN and model.modelFile:
            file_format = None
            algorithm = model.trainingStrategy.algorithm
            for format_ in algorithm.fileFormats.all():
                if validated_data['file'].name.endswith(format_.fileExtension):
                    file_format = format_
                    break
            if not file_format:
                raise self.InvalidFileFormatError(f"The extension for file '{validated_data['file'].name}' of the submitted file did not match any of the known formats for algorithm: ({algorithm.name}).")

            if model.modelFile.format.fileExtension == file_format.fileExtension:
                model.modelFile.file.save(model.modelFile.path, validated_data['file'])
            else:
                model.modelFile.delete()
                ModelFile.objects.create(
                    model=model,
                    kind=ModelFile.MAIN,
                    format=file_format
                )

            return model.modelFile
        else:
            ModelFileFormat = modelling.models.ModelFileFormat
            file_format = None
            for format_ in ModelFileFormat.objects.all():
                if validated_data['file'].name.endswith(format_.fileExtension):
                    file_format = format_
                    break
            ret = ModelFile.objects.create(
                modelInstance=model,
                kind=validated_data['kind'],
                format=file_format if file_format else ModelFileFormat.objects.get_or_create(
                    fileExtension=validated_data['file'].name.split('.')[-1]
                )
            )
            if ret.kind == ModelFile.MAIN:
                ret.file.save(ret.generateMainFileName(model, file_format), validated_data['file'])
            else:
                ret.file.save(ret.generateAuxFileName(model, validated_data['file'].name), validated_data['file'])

            return ret

class ModelSerializer(serializers.HyperlinkedModelSerializer):
    project = serializers.PrimaryKeyRelatedField(many=False, queryset=Project.objects.all())
    performance = ModelPerformanceSerializer(many=True, read_only=True)
    trainingStrategy = TrainingStrategySerializer(many=False)
    validationStrategy = BasicValidationStrategyInitSerializer(many=False, required=False)
    build = serializers.BooleanField(default=True)
    taskID = serializers.UUIDField(required=False, read_only=True, allow_null=True)
    modelFile = ModelFileSerializer(many=False, read_only=True, allow_null=True)

    def __init__(self, *args, builder_class=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.builder_class = self.instance.builder.name if self.instance and isinstance(self.instance, modelling.models.Model) else builder_class

    # def validate(self, attrs):
    #     if "modelFile" in attrs and "validationStrategy" in attrs:
    #         raise ValidationError("If 'modelFile' is present, 'validationStrategy' field should be empty.")
    #     if "modelFile" not in attrs and "validationStrategy" not in attrs:
    #         raise ValidationError("You have to specify 'modelFile' if you omit 'validationStrategy'.")
    #     return super().validate(attrs)

    @staticmethod
    def saveParameters(strat_instance, strat_data):
        if 'parameters' not in strat_data:
            return
        for param_name in strat_data['parameters']:
            parameter = modelling.models.ModelParameter.objects.get(
                name=param_name
                , algorithm__name=strat_data['algorithm'].name
            )
            value_class = modelling.models.PARAM_VALUE_CTYPE_TO_MODEL_MAP[parameter.contentType]
            parameter_value = value_class(
                parameter=parameter
                , strategy=strat_instance
                , value=value_class.parseValue(strat_data['parameters'][param_name]))
            parameter_value.save()

    class Meta:
        model = modelling.models.Model
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'trainingStrategy', 'validationStrategy', 'performance', 'modelFile', 'build', 'taskID')
        read_only_fields = ('id', 'created', 'updated', 'performance', 'modelFile', 'taskID')

    def useBuilder(self, builder_class):
        self.builder_class = builder_class

    def create(self, validated_data, **kwargs):
        instance = self.Meta.model.objects.create(
                name=validated_data['name'],
                description=validated_data['description'],
                project=validated_data['project'],
                builder=modelling.models.ModelBuilder.objects.get_or_create(
                    name=self.builder_class if type(self.builder_class) == str else self.builder_class.__name__
                )[0],
                **kwargs
        )
        instance.build = validated_data["build"]
        return instance

class BasicValidationStrategySerializer(BasicValidationStrategyInitSerializer):
    metrics = ModelPerformanceMetricSerializer(many=True)