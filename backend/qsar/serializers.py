"""
serializers

Created by: Martin Sicho
On: 13-01-20, 11:07
"""

from rest_framework import serializers

from commons.serializers import GenericModelSerializerMixIn
from projects.models import Project
from . import models

class ModelPerformanceMetricSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = models.ModelPerformanceMetric
        fields = ('name', 'description')

class ModelPerformanceSerializer(GenericModelSerializerMixIn, serializers.HyperlinkedModelSerializer):
    className = GenericModelSerializerMixIn.className
    extraArgs = GenericModelSerializerMixIn.extraArgs
    metric = ModelPerformanceMetricSerializer(many=False)

    class Meta:
        model = models.ModelPerformance
        fields = ('id', 'value', 'metric', 'className', 'extraArgs')

class ModelFileFormatSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = models.ModelFileFormat
        fields = ('fileExtension', 'description')

class ModelParameterSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = models.ModelParameter
        fields = ('name', 'contentType')

class AlgorithmModeSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = models.AlgorithmMode
        fields = ('name',)

class AlgorithmSerializer(serializers.HyperlinkedModelSerializer):
    fileFormats = ModelFileFormatSerializer(many=True)
    parameters = ModelParameterSerializer(many=True)
    validModes = AlgorithmModeSerializer(many=True)

    class Meta:
        model = models.Algorithm
        fields = ('id', 'name', 'fileFormats', 'parameters', 'validModes')

class ModelParameterValueSerializer(serializers.HyperlinkedModelSerializer):
    parameter = ModelParameterSerializer(many=False)
    value = serializers.CharField()

    class Meta:
        model = models.ModelParameterValue
        fields = ('parameter', 'value')

class TrainingStrategySerializer(serializers.HyperlinkedModelSerializer):
    algorithm = AlgorithmSerializer(many=False)
    parameters = ModelParameterValueSerializer(many=True)
    mode = AlgorithmModeSerializer(many=False)

    class Meta:
        model = models.TrainingStrategy
        fields = ('algorithm', 'mode', 'parameters')

class DescriptorGroupSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = models.DescriptorGroup
        fields = ('name',)

class QSARTrainingStrategySerializer(TrainingStrategySerializer):
    descriptors = DescriptorGroupSerializer(many=True)

    class Meta:
        model = models.QSARTrainingStrategy
        fields = TrainingStrategySerializer.Meta.fields + ('descriptors', 'activityThreshold')

class QSARTrainingStrategySerializerInit(QSARTrainingStrategySerializer):
    descriptors = serializers.PrimaryKeyRelatedField(many=True, queryset=models.DescriptorGroup.objects.all())
    algorithm = serializers.PrimaryKeyRelatedField(many=False, queryset=models.Algorithm.objects.all())
    parameters = serializers.DictField(allow_empty=True, child=serializers.CharField())
    mode = serializers.PrimaryKeyRelatedField(many=False, queryset=models.AlgorithmMode.objects.all())

    class Meta:
        model = models.QSARTrainingStrategy
        fields = QSARTrainingStrategySerializer.Meta.fields

class BasicValidationStrategySerializer(serializers.HyperlinkedModelSerializer):
    metrics = serializers.PrimaryKeyRelatedField(many=True, queryset=models.ModelPerformanceMetric.objects.all())
    cvFolds = serializers.IntegerField(min_value=0)
    validSetSize = serializers.FloatField(min_value=0)

    class Meta:
        model = models.BasicValidationStrategy
        fields = ('cvFolds', 'validSetSize', 'metrics')

class ModelSerializer(serializers.HyperlinkedModelSerializer):
    project = serializers.PrimaryKeyRelatedField(many=False, queryset=Project.objects.all())
    performance = ModelPerformanceSerializer(many=True)
    trainingStrategy = TrainingStrategySerializer(many=False)
    validationStrategy = BasicValidationStrategySerializer(many=False)

    class Meta:
        model = models.Model
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'trainingStrategy', 'validationStrategy', 'performance')
        read_only_fields = ('id', 'created', 'updated', 'performance')

class QSARModelSerializer(ModelSerializer):
    trainingStrategy = QSARTrainingStrategySerializer(many=False)
    molset = serializers.PrimaryKeyRelatedField(many=False, queryset=models.MolSet.objects.all())
    predictions = serializers.PrimaryKeyRelatedField(many=True, queryset=models.ActivitySet.objects.all())
    taskID = serializers.UUIDField(required=False)

    class Meta:
        model = models.QSARModel
        fields = ModelSerializer.Meta.fields + ('molset', 'predictions', 'taskID')
        read_only_fields = ModelSerializer.Meta.read_only_fields + ('predictions', 'taskID')

class QSARModelSerializerInit(QSARModelSerializer):
    trainingStrategy = QSARTrainingStrategySerializerInit(many=False)

    class Meta:
        model = models.QSARModel
        fields = ('name', 'description', 'project', 'trainingStrategy', 'validationStrategy', 'molset')
        read_only_fields = None

    def create(self, validated_data):
        instance = models.QSARModel.objects.create(
            name=validated_data['name'],
            description=validated_data['description'],
            project=validated_data['project'],
            molset=validated_data['molset'],
        )
        models.ModelActivitySet.objects.create(model=instance, project=instance.project)

        strat_data = validated_data['trainingStrategy']
        trainingStrategy = models.QSARTrainingStrategy(
            modelInstance = instance,
            algorithm = strat_data['algorithm'],
            mode = strat_data['mode'],
            activityThreshold = strat_data['activityThreshold']
        )
        trainingStrategy.save()
        trainingStrategy.descriptors.set(strat_data['descriptors'])
        trainingStrategy.save()

        for param_name in strat_data['parameters']:
            parameter = models.ModelParameter.objects.get(
                name=param_name
                , algorithm__name=strat_data['algorithm'].name
            )
            value_class = models.PARAM_VALUE_CTYPE_TO_MODEL_MAP[parameter.contentType]
            parameter_value = value_class(
                parameter=parameter
                , strategy=trainingStrategy
                , value=value_class.parseValue(strat_data['parameters'][param_name]))
            parameter_value.save()

        strat_data = validated_data['validationStrategy']
        validationStrategy = models.BasicValidationStrategy.objects.create(
            modelInstance = instance,
            cvFolds=strat_data['cvFolds'],
            validSetSize=strat_data['validSetSize']
        )
        validationStrategy.metrics.set(strat_data['metrics'])
        validationStrategy.save()

        return instance

