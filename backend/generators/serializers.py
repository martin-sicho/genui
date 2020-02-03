"""
serializers

Created by: Martin Sicho
On: 27-01-20, 17:00
"""
from rest_framework import serializers

from compounds.models import MolSet
from compounds.serializers import MolSetSerializer
from modelling.models import ModelPerformanceMetric
from modelling.serializers import ModelSerializer, ValidationStrategySerializer, TrainingStrategySerializer, \
    TrainingStrategyInitSerializer
from qsar.models import QSARModel
from qsar.serializers import QSARModelSerializer
from . import models

class DrugExCorpusSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = models.DrugeExCorpus
        fields = ("corpusFile", "vocFile")
        read_only_fields = ("corpusFile", "vocFile")

class DrugExValidationStrategySerializer(ValidationStrategySerializer):

    class Meta:
        model = models.DrugExValidationStrategy
        fields = ValidationStrategySerializer.Meta.fields + ("validSetSize",)

class DrugExValidationStrategyInitSerializer(DrugExValidationStrategySerializer):
    metrics = serializers.PrimaryKeyRelatedField(many=True, queryset=ModelPerformanceMetric.objects.all())

    class Meta:
        model = models.DrugExValidationStrategy
        fields = DrugExValidationStrategySerializer.Meta.fields + ("validSetSize",)

class DrugExTrainingStrategySerializer(TrainingStrategySerializer):

    class Meta:
        model = models.DrugExNetTrainingStrategy
        fields = TrainingStrategySerializer.Meta.fields

class DrugExTrainingStrategyInitSerializer(TrainingStrategyInitSerializer):

    class Meta:
        model = models.DrugExNetTrainingStrategy
        fields = TrainingStrategyInitSerializer.Meta.fields


class DrugExNetSerializer(ModelSerializer):
    molset = MolSetSerializer(many=False)
    corpus = DrugExCorpusSerializer(many=False, read_only=True)
    trainingStrategy = DrugExTrainingStrategySerializer(many=False)
    validationStrategy = DrugExValidationStrategySerializer(many=False)
    parent = serializers.PrimaryKeyRelatedField(many=False, queryset=models.DrugExNet.objects.all(), required=False)

    class Meta:
        model = models.DrugExNet
        fields = [x for x in ModelSerializer.Meta.fields if x not in ('performance',)] + ["molset", "corpus", "parent"]
        read_only_fields = [x for x in ModelSerializer.Meta.read_only_fields if x not in ('performance',)] + ["corpus"]


class DrugExNetInitSerializer(DrugExNetSerializer):
    molset = serializers.PrimaryKeyRelatedField(many=False, queryset=MolSet.objects.all(), required=True)
    trainingStrategy = DrugExTrainingStrategyInitSerializer(many=False)
    validationStrategy = DrugExValidationStrategyInitSerializer(many=False)

    class Meta:
        model = models.DrugExNet
        fields = DrugExNetSerializer.Meta.fields
        read_only_fields = DrugExNetSerializer.Meta.read_only_fields

    def create(self, validated_data):
        instance = models.DrugExNet.objects.create(
            name=validated_data['name'],
            description=validated_data['description'],
            project=validated_data['project'],
            molset=validated_data['molset'],
        )
        if "parent" in validated_data:
            instance.parent=models.DrugExNet.objects.get(pk=validated_data['parent'])
            instance.save()

        strat_data = validated_data['trainingStrategy']
        trainingStrategy = models.DrugExNetTrainingStrategy.objects.create(
            modelInstance=instance,
            algorithm=strat_data['algorithm'],
            mode=strat_data['mode'],
        )

        self.saveParameters(trainingStrategy, strat_data)

        strat_data = validated_data['validationStrategy']
        validationStrategy = models.DrugExValidationStrategy.objects.create(
            modelInstance = instance,
            validSetSize=strat_data['validSetSize']
        )
        strat_data.update({
            'metrics' : [x for x in ModelPerformanceMetric.objects.filter(name__in=('SMILES_ER', 'DrExLoss'))] # TODO: every validation strategy should have a set of acceptable validation metrics so that this doesn't have to be hardcoded
        })
        validationStrategy.metrics.set(strat_data['metrics'])
        validationStrategy.save()

        return instance

class DrugExAgentTrainingStrategySerializer(TrainingStrategySerializer):

    class Meta:
        model = models.DrugExAgentTrainingStrategy
        fields = TrainingStrategySerializer.Meta.fields

class DrugExAgentTrainingStrategyInitSerializer(TrainingStrategyInitSerializer):

    class Meta:
        model = models.DrugExAgentTrainingStrategy
        fields = TrainingStrategyInitSerializer.Meta.fields

class DrugExAgentValidationStrategySerializer(ValidationStrategySerializer):

    class Meta:
        model = models.DrugExAgentValidationStrategy
        fields = ValidationStrategySerializer.Meta.fields

class DrugExAgentValidationStrategyInitSerializer(DrugExAgentValidationStrategySerializer):

    class Meta:
        model = models.DrugExAgentValidationStrategy
        fields = ValidationStrategySerializer.Meta.fields

class DrugExAgentSerializer(ModelSerializer):
    trainingStrategy = DrugExAgentTrainingStrategySerializer(many=False)
    validationStrategy = DrugExAgentValidationStrategySerializer(many=False)
    environment = QSARModelSerializer(many=False)
    explorationNet = DrugExNetSerializer(many=False)
    exploitationNet = DrugExNetSerializer(many=False)

    class Meta:
        model = models.DrugExAgent
        fields = [x for x in ModelSerializer.Meta.fields if x not in ('performance',)]  + ["environment", "explorationNet", "exploitationNet"]
        read_only_fields = [x for x in ModelSerializer.Meta.read_only_fields if x not in ('performance',)]

class DrugExAgentInitSerializer(DrugExAgentSerializer):
    trainingStrategy = DrugExAgentTrainingStrategyInitSerializer(many=False)
    validationStrategy = DrugExAgentValidationStrategyInitSerializer(many=False, required=False)
    environment = serializers.PrimaryKeyRelatedField(many=False, queryset=QSARModel.objects.all())
    explorationNet = serializers.PrimaryKeyRelatedField(many=False, queryset=models.DrugExNet.objects.all())
    exploitationNet = serializers.PrimaryKeyRelatedField(many=False, queryset=models.DrugExNet.objects.all())

    class Meta:
        model = models.DrugExAgent
        fields = DrugExAgentSerializer.Meta.fields
        read_only_fields = DrugExAgentSerializer.Meta.read_only_fields

    def create(self, validated_data):
        instance = models.DrugExAgent.objects.create(
            name=validated_data['name'],
            description=validated_data['description'],
            project=validated_data['project'],
            environment=validated_data['environment'],
            exploitationNet=validated_data['exploitationNet'],
            explorationNet=validated_data['explorationNet'],
        )

        strat_data = validated_data['trainingStrategy']
        trainingStrategy = models.DrugExAgentTrainingStrategy.objects.create(
            modelInstance=instance,
            algorithm=strat_data['algorithm'],
            mode=strat_data['mode'],
        )
        self.saveParameters(trainingStrategy, strat_data)

        validationStrategy = models.DrugExValidationStrategy.objects.create(
                modelInstance = instance
            )
        if 'validationStrategy' in validated_data:
            strat_data = validated_data['validationStrategy']
        else:
            strat_data = {
                'metrics' : [x for x in ModelPerformanceMetric.objects.filter(name__in=('SMILES_ER', 'SMILES_UQR', 'DrExActivity'))] # TODO: every validation strategy should have a set of acceptable validation metrics so that this doesn't have to be hardcoded
            }
        validationStrategy.metrics.set(strat_data['metrics'])
        validationStrategy.save()

        return instance
