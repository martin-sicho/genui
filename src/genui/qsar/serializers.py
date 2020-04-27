"""
serializers

Created by: Martin Sicho
On: 13-01-20, 11:07
"""

from rest_framework import serializers

from genui.compounds.models import ActivityTypes, ActivitySet, ActivityUnits
from genui.compounds.serializers import MolSetSerializer, ActivitySetSerializer, ActivityTypeSerializer, \
    ActivityUnitsSerializer
from genui.modelling.serializers import TrainingStrategySerializer, BasicValidationStrategyInitSerializer, ModelSerializer, \
    BasicValidationStrategySerializer, TrainingStrategyInitSerializer
from . import models
from genui.modelling.models import BasicValidationStrategy


class DescriptorGroupSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = models.DescriptorGroup
        fields = ('id', 'name',)

class QSARTrainingStrategySerializer(TrainingStrategySerializer):
    descriptors = DescriptorGroupSerializer(many=True)
    activityType = ActivityTypeSerializer(many=False)
    activitySet = ActivitySetSerializer(many=False)

    class Meta:
        model = models.QSARTrainingStrategy
        fields = TrainingStrategySerializer.Meta.fields + ('descriptors', 'activityThreshold', 'activityType', 'activitySet')

class QSARTrainingStrategyInitSerializer(TrainingStrategyInitSerializer):
    descriptors = serializers.PrimaryKeyRelatedField(many=True, queryset=models.DescriptorGroup.objects.all(), allow_empty=False)
    activityType = serializers.PrimaryKeyRelatedField(many=False, queryset=ActivityTypes.objects.all(), required=False)
    activitySet = serializers.PrimaryKeyRelatedField(many=False, queryset=ActivitySet.objects.all(), required=False)

    class Meta:
        model = models.QSARTrainingStrategy
        fields = QSARTrainingStrategySerializer.Meta.fields


class QSARModelSerializer(ModelSerializer):
    trainingStrategy = QSARTrainingStrategySerializer(many=False)
    validationStrategy = BasicValidationStrategySerializer(many=False, required=False)
    molset = MolSetSerializer(many=False, required=False)
    predictions = serializers.PrimaryKeyRelatedField(many=True, queryset=models.ActivitySet.objects.all())
    predictionsType = ActivityTypeSerializer(many=False)
    predictionsUnits = ActivityUnitsSerializer(many=False)

    class Meta:
        model = models.QSARModel
        fields = ModelSerializer.Meta.fields + ('molset', 'predictions', 'predictionsType', 'predictionsUnits')
        read_only_fields = ModelSerializer.Meta.read_only_fields + ('predictions',)

class QSARModelInitSerializer(QSARModelSerializer):
    trainingStrategy = QSARTrainingStrategyInitSerializer(many=False)
    validationStrategy = BasicValidationStrategyInitSerializer(many=False, required=False)
    molset = serializers.PrimaryKeyRelatedField(many=False, queryset=models.MolSet.objects.all(), required=False)
    predictionsType = serializers.CharField(required=False, max_length=128, allow_null=False)
    predictionsUnits = serializers.CharField(required=False, max_length=128, allow_null=True)

    class Meta:
        model = models.QSARModel
        fields = [x for x in QSARModelSerializer.Meta.fields if x not in ('predictions',)]
        read_only_fields = QSARModelSerializer.Meta.read_only_fields

    def is_valid(self, raise_exception=True):
        ret = super().is_valid(raise_exception)
        data = self.validated_data
        tr_strat_data = data['trainingStrategy']

        if data['build'] and tr_strat_data['mode'].name == "classification" and ('activityThreshold' not in  tr_strat_data or tr_strat_data['activityThreshold'] is None):
            raise serializers.ValidationError("You must specify an activity threshold for a classification model.")

        if data['build'] and ('activityType' not in  tr_strat_data or tr_strat_data['activityType'] is None):
            raise serializers.ValidationError("You have to specify the activity type of the training data. Use the 'activityType' parameter in 'trainingStrategy'.")

        if data['build'] and ('activitySet' not in  tr_strat_data or tr_strat_data['activitySet'] is None):
            raise serializers.ValidationError("You have to specify the activity set that contains the true activities for training. Use the 'activitySet' parameter in 'trainingStrategy'.")

        if not data["build"] and ("predictionsType" not in data or "predictionsUnits" not in data or not data["predictionsType"]):
            raise serializers.ValidationError("You have to specify the type and units of the predicted values if you are not building the model from existing data. Both 'predictionsType' and 'predictionsUnits' must be specified. You can set 'predictionsUnits' to 'null' if the model output variable has no dimension.")

        return ret

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
            activityThreshold = strat_data['activityThreshold'] if 'activityThreshold' in strat_data else None,
            activitySet=strat_data['activitySet'] if 'activitySet' in strat_data else None,
            activityType=strat_data['activityType'] if 'activityType' in strat_data else None
        )
        trainingStrategy.save()
        trainingStrategy.descriptors.set(strat_data['descriptors'])
        trainingStrategy.save()

        self.saveParameters(trainingStrategy, strat_data)

        if 'validationStrategy' in validated_data:
            strat_data = validated_data['validationStrategy']
            validationStrategy = BasicValidationStrategy.objects.create(
                modelInstance = instance,
                cvFolds=strat_data['cvFolds'],
                validSetSize=strat_data['validSetSize']
            )
            validationStrategy.metrics.set(strat_data['metrics'])
            validationStrategy.save()

        if "predictionsType" in validated_data:
            instance.predictionsType = ActivityTypes.objects.get_or_create(
                value=validated_data["predictionsType"]
            )[0]

        if "predictionsUnits" in validated_data and validated_data["predictionsUnits"]:
            instance.predictionsUnits = ActivityUnits.objects.get_or_create(
                value=validated_data["predictionsUnits"]
            )[0]

        instance.save()

        return instance

class ModelActivitySetSerializer(ActivitySetSerializer):
    model = serializers.PrimaryKeyRelatedField(many=False, queryset=models.QSARModel.objects.all())
    taskID = serializers.UUIDField(read_only=True, required=False)

    class Meta:
        model = models.ModelActivitySet
        fields = ActivitySetSerializer.Meta.fields + ('model', 'taskID')
        read_only_fields = ActivitySetSerializer.Meta.read_only_fields + ('taskID', 'model', 'project')