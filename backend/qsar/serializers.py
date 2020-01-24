"""
serializers

Created by: Martin Sicho
On: 13-01-20, 11:07
"""

from rest_framework import serializers

import modelling.models
from compounds.serializers import MolSetSerializer
from modelling.serializers import TrainingStrategySerializer, BasicValidationStrategySerializer, ModelSerializer, \
    ValidationStrategySerializer
from . import models


class DescriptorGroupSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = models.DescriptorGroup
        fields = ('id', 'name',)

class QSARTrainingStrategySerializer(TrainingStrategySerializer):
    descriptors = DescriptorGroupSerializer(many=True)

    class Meta:
        model = models.QSARTrainingStrategy
        fields = TrainingStrategySerializer.Meta.fields + ('descriptors', 'activityThreshold')

class QSARTrainingStrategyInitSerializer(QSARTrainingStrategySerializer):
    descriptors = serializers.PrimaryKeyRelatedField(many=True, queryset=models.DescriptorGroup.objects.all())
    algorithm = serializers.PrimaryKeyRelatedField(many=False, queryset=modelling.models.Algorithm.objects.all())
    parameters = serializers.DictField(allow_empty=True, child=serializers.CharField())
    mode = serializers.PrimaryKeyRelatedField(many=False, queryset=modelling.models.AlgorithmMode.objects.all())

    class Meta:
        model = models.QSARTrainingStrategy
        fields = QSARTrainingStrategySerializer.Meta.fields


class QSARModelSerializer(ModelSerializer):
    trainingStrategy = QSARTrainingStrategySerializer(many=False)
    validationStrategy = ValidationStrategySerializer(many=False)
    molset = MolSetSerializer(many=False)
    predictions = serializers.PrimaryKeyRelatedField(many=True, queryset=models.ActivitySet.objects.all())
    taskID = serializers.UUIDField(required=False)

    class Meta:
        model = models.QSARModel
        fields = ModelSerializer.Meta.fields + ('molset', 'predictions', 'taskID', 'modelFile')
        read_only_fields = ModelSerializer.Meta.read_only_fields + ('predictions', 'taskID', 'modelFile')

class QSARModelInitSerializer(QSARModelSerializer):
    trainingStrategy = QSARTrainingStrategyInitSerializer(many=False)
    validationStrategy = BasicValidationStrategySerializer(many=False)
    molset = serializers.PrimaryKeyRelatedField(many=False, queryset=models.MolSet.objects.all())

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
            parameter = modelling.models.ModelParameter.objects.get(
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
        validationStrategy = modelling.models.BasicValidationStrategy.objects.create(
            modelInstance = instance,
            cvFolds=strat_data['cvFolds'],
            validSetSize=strat_data['validSetSize']
        )
        validationStrategy.metrics.set(strat_data['metrics'])
        validationStrategy.save()

        return instance

