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
    parameterValues = ModelParameterValueSerializer(many=True)
    metrics = ModelPerformanceMetricSerializer(many=True)
    mode = AlgorithmModeSerializer(many=False)

    class Meta:
        model = models.TrainingStrategy
        fields = ('algorithm', 'parameterValues', 'metrics', 'mode')

class DescriptorGroupSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = models.DescriptorGroup
        fields = ('name',)

class QSARTrainingStrategySerializer(TrainingStrategySerializer):
    descriptors = DescriptorGroupSerializer(many=True)

    class Meta:
        model = models.QSARTrainingStrategy
        fields = TrainingStrategySerializer.Meta.fields + ('descriptors', 'activityThreshold')

class ValidationStrategySerializer(GenericModelSerializerMixIn, serializers.HyperlinkedModelSerializer):
    className = GenericModelSerializerMixIn.className
    extraArgs = GenericModelSerializerMixIn.extraArgs

    class Meta:
        model = models.ValidationStrategy
        fields = ('className', 'extraArgs')

class ModelSerializer(serializers.HyperlinkedModelSerializer):
    project = serializers.PrimaryKeyRelatedField(many=False, queryset=Project.objects.all())
    performance = ModelPerformanceSerializer(many=True)
    trainingStrategy = TrainingStrategySerializer(many=False)
    validationStrategy = ValidationStrategySerializer(many=False)

    class Meta:
        model = models.Model
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'trainingStrategy', 'validationStrategy', 'performance')
        read_only_fields = ('id', 'created', 'updated', 'performance')

class QSARModelSerializer(ModelSerializer):
    trainingStrategy = QSARTrainingStrategySerializer(many=False)
    molset = serializers.PrimaryKeyRelatedField(many=False, queryset=models.MolSet.objects.all())

    class Meta:
        model = models.QSARModel
        fields = ModelSerializer.Meta.fields + ('molset', 'activities')
        read_only_fields = ModelSerializer.Meta.read_only_fields + ('activities',)

