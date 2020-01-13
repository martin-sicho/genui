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

class AlgorithmSerializer(serializers.HyperlinkedModelSerializer):
    fileFormats = ModelFileFormatSerializer(many=True)
    parameters = ModelParameterSerializer(many=True)

    class Meta:
        model = models.Algorithm
        fields = ('id', 'name', 'fileFormats', 'parameters')

class ModelParameterValueSerializer(GenericModelSerializerMixIn, serializers.HyperlinkedModelSerializer):
    parameter = ModelParameterSerializer(many=False)
    className = GenericModelSerializerMixIn.className
    extraArgs = GenericModelSerializerMixIn.extraArgs

    class Meta:
        model = models.ModelParameterValue
        fields = ('id', 'parameter', 'className', 'extraArgs')

class TrainingStrategySerializer(GenericModelSerializerMixIn, serializers.HyperlinkedModelSerializer):
    className = GenericModelSerializerMixIn.className
    extraArgs = GenericModelSerializerMixIn.extraArgs
    algorithm = AlgorithmSerializer(many=False)
    parameters = ModelParameterValueSerializer(many=True)
    fileFormat = ModelFileFormatSerializer(many=False)
    metrics = ModelPerformanceMetricSerializer(many=True)

    class Meta:
        model = models.TrainingStrategy
        fields = ('id', 'algorithm', 'parameters', 'fileFormat', 'metrics', 'className', 'extraArgs')

class ValidationStrategySerializer(GenericModelSerializerMixIn, serializers.HyperlinkedModelSerializer):
    className = GenericModelSerializerMixIn.className
    extraArgs = GenericModelSerializerMixIn.extraArgs

    class Meta:
        model = models.ValidationStrategy
        fields = ('id', 'className', 'extraArgs')

class ModelSerializer(serializers.HyperlinkedModelSerializer):
    project = serializers.PrimaryKeyRelatedField(many=False, queryset=Project.objects.all())
    performance = ModelPerformanceSerializer(many=True)
    trainingStrategy = TrainingStrategySerializer(many=False)
    validationStrategy = ValidationStrategySerializer(many=False)

    class Meta:
        model = models.Model
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'trainingStrategy', 'validationStrategy', 'performance')
        read_only_fields = ('id', 'created', 'updated', 'performance')

class DescriptorGroupSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = models.DescriptorGroup
        fields = ('id', 'name')

class QSARModelSerializer(ModelSerializer):
    descriptors = DescriptorGroupSerializer(many=True)

    class Meta:
        model = models.QSARModel
        fields = ModelSerializer.Meta.fields + ('molset', 'activities', 'descriptors')
        read_only_fields = ModelSerializer.Meta.read_only_fields + ('activities',)

