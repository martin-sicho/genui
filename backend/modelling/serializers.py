"""
serializers

Created by: Martin Sicho
On: 24-01-20, 14:44
"""
from rest_framework import serializers

from commons.serializers import GenericModelSerializerMixIn
from projects.models import Project
from qsar import models


class ModelPerformanceMetricSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = models.ModelPerformanceMetric
        fields = ('id', 'name', 'description')


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
        fields = ('id', 'fileExtension', 'description')


class ModelParameterSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = models.ModelParameter
        fields = ('id', 'name', 'contentType')


class AlgorithmModeSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = models.AlgorithmMode
        fields = ('id', 'name',)


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
        fields = ('id','parameter', 'value')


class TrainingStrategySerializer(serializers.HyperlinkedModelSerializer):
    algorithm = AlgorithmSerializer(many=False)
    parameters = ModelParameterValueSerializer(many=True)
    mode = AlgorithmModeSerializer(many=False)

    class Meta:
        model = models.TrainingStrategy
        fields = ('algorithm', 'mode', 'parameters')


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


class ValidationStrategySerializer(BasicValidationStrategySerializer):
    metrics = ModelPerformanceMetricSerializer(many=True)