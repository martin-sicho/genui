"""
serializers

Created by: Martin Sicho
On: 13-01-20, 11:07
"""

from rest_framework import serializers
from rest_framework.exceptions import ValidationError

import modelling.models
from compounds.serializers import MolSetSerializer, ActivitySetSerializer
from modelling.serializers import TrainingStrategySerializer, BasicValidationStrategyInitSerializer, ModelSerializer, \
    BasicValidationStrategySerializer, TrainingStrategyInitSerializer
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

class QSARTrainingStrategyInitSerializer(TrainingStrategyInitSerializer):
    descriptors = serializers.PrimaryKeyRelatedField(many=True, queryset=models.DescriptorGroup.objects.all(), allow_empty=False)

    class Meta:
        model = models.QSARTrainingStrategy
        fields = QSARTrainingStrategySerializer.Meta.fields


class QSARModelSerializer(ModelSerializer):
    trainingStrategy = QSARTrainingStrategySerializer(many=False)
    validationStrategy = BasicValidationStrategySerializer(many=False, required=False)
    molset = MolSetSerializer(many=False, required=False)
    predictions = serializers.PrimaryKeyRelatedField(many=True, queryset=models.ActivitySet.objects.all())

    class Meta:
        model = models.QSARModel
        fields = ModelSerializer.Meta.fields + ('molset', 'predictions')
        read_only_fields = ModelSerializer.Meta.read_only_fields + ('predictions',)

    def validate(self, attrs):
        attrs = super().validate(attrs)
        # if 'modelFile' not in attrs and 'molset' not in attrs:
        #     raise ValidationError("You must specify 'modelFile' if there is no training set.")
        if 'molset' in attrs and 'activityThreshold' not in attrs['trainingStrategy']:
            raise ValidationError("You have to specify 'activityThreshold' if you also specify 'molset'.")
        return attrs


class QSARModelInitSerializer(QSARModelSerializer):
    trainingStrategy = QSARTrainingStrategyInitSerializer(many=False)
    validationStrategy = BasicValidationStrategyInitSerializer(many=False, required=False)
    molset = serializers.PrimaryKeyRelatedField(many=False, queryset=models.MolSet.objects.all(), required=False)

    class Meta:
        model = models.QSARModel
        fields = ModelSerializer.Meta.fields + ('molset',)
        read_only_fields = ModelSerializer.Meta.read_only_fields

    def create(self, validated_data, **kwargs):
        instance = super().create(
            validated_data
            , molset=validated_data['molset'] if 'molset' in validated_data else None
            , **kwargs
        )

        strat_data = validated_data['trainingStrategy']
        trainingStrategy = models.QSARTrainingStrategy(
            modelInstance = instance,
            algorithm = strat_data['algorithm'],
            mode = strat_data['mode'],
            activityThreshold = strat_data['activityThreshold'] if 'activityThreshold' in strat_data else None
        )
        trainingStrategy.save()
        trainingStrategy.descriptors.set(strat_data['descriptors'])
        trainingStrategy.save()

        self.saveParameters(trainingStrategy, strat_data)

        if 'validationStrategy' in validated_data:
            strat_data = validated_data['validationStrategy']
            validationStrategy = modelling.models.BasicValidationStrategy.objects.create(
                modelInstance = instance,
                cvFolds=strat_data['cvFolds'],
                validSetSize=strat_data['validSetSize']
            )
            validationStrategy.metrics.set(strat_data['metrics'])
            validationStrategy.save()

        return instance

class ModelActivitySetSerializer(ActivitySetSerializer):
    model = serializers.PrimaryKeyRelatedField(many=False, queryset=models.QSARModel.objects.all())
    taskID = serializers.UUIDField(read_only=True, required=False)

    class Meta:
        model = models.ModelActivitySet
        fields = ActivitySetSerializer.Meta.fields + ('model', 'taskID')
        read_only_fields = ActivitySetSerializer.Meta.read_only_fields + ('taskID', 'model', 'project')