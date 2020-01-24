"""
serializers

Created by: Martin Sicho
On: 24-01-20, 14:44
"""
from rest_framework import serializers

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


class BasicValidationStrategySerializer(serializers.HyperlinkedModelSerializer):
    metrics = serializers.PrimaryKeyRelatedField(many=True, queryset=modelling.models.ModelPerformanceMetric.objects.all())
    cvFolds = serializers.IntegerField(min_value=0)
    validSetSize = serializers.FloatField(min_value=0)

    class Meta:
        model = modelling.models.BasicValidationStrategy
        fields = ('cvFolds', 'validSetSize', 'metrics')


class ModelSerializer(serializers.HyperlinkedModelSerializer):
    project = serializers.PrimaryKeyRelatedField(many=False, queryset=Project.objects.all())
    performance = ModelPerformanceSerializer(many=True)
    trainingStrategy = TrainingStrategySerializer(many=False)
    validationStrategy = BasicValidationStrategySerializer(many=False)

    class Meta:
        model = modelling.models.Model
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'trainingStrategy', 'validationStrategy', 'performance')
        read_only_fields = ('id', 'created', 'updated', 'performance')


class ValidationStrategySerializer(BasicValidationStrategySerializer):
    metrics = ModelPerformanceMetricSerializer(many=True)