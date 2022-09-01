"""
serializers

Created by: Martin Sicho
On: 5/3/20, 6:43 PM
"""
from django.db.models import Q
from rest_framework import serializers

from genui.projects.models import Project
from genui.qsar.models import QSARModel
from . import models
from genui.compounds.models import MolSet
from genui.compounds.serializers import MolSetSerializer
from genui.generators.serializers import GeneratorSerializer
from genui.models.genuimodels.bases import Algorithm
from genui.models.models import ModelPerformanceMetric, Model
from genui.models.serializers import ValidationStrategySerializer, TrainingStrategySerializer, \
    TrainingStrategyInitSerializer, ModelSerializer, ValidationStrategyInitSerializer
# from genui.qsar.models import QSARModel
# from genui.qsar.serializers import QSARModelSerializer


class DrugExValidationStrategySerializer(ValidationStrategySerializer):

    class Meta:
        model = models.DrugExNetValidation
        fields = ValidationStrategySerializer.Meta.fields + ("validSetSize",)


class DrugExValidationStrategyInitSerializer(ValidationStrategyInitSerializer):
    metrics = serializers.PrimaryKeyRelatedField(many=True, queryset=ModelPerformanceMetric.objects.all(), required=False)

    class Meta:
        model = models.DrugExNetValidation
        fields = DrugExValidationStrategySerializer.Meta.fields + ("validSetSize",)


class DrugExTrainingStrategySerializer(TrainingStrategySerializer):

    class Meta:
        model = models.DrugExNetTraining
        fields = TrainingStrategySerializer.Meta.fields + ("inputType", "modelClass")


class DrugExTrainingStrategyInitSerializer(TrainingStrategyInitSerializer):

    class Meta:
        model = models.DrugExNetTraining
        fields = TrainingStrategyInitSerializer.Meta.fields + ("inputType", "modelClass")


class DrugExNetSerializer(ModelSerializer):
    molset = MolSetSerializer(many=False, required=False)
    trainingStrategy = DrugExTrainingStrategySerializer(many=False)
    validationStrategy = DrugExValidationStrategySerializer(many=False, required=False)
    parent = serializers.SerializerMethodField("get_parent")

    class Meta:
        model = models.DrugExNet
        fields = [x for x in ModelSerializer.Meta.fields if x not in ('performance',)] + ["molset", "parent"]
        read_only_fields = [x for x in ModelSerializer.Meta.read_only_fields if x not in ('performance',)]

    def get_parent(self, obj):
        if obj.parent is not None:
            return DrugExNetSerializer(obj.parent).data
        else:
            return None


class DrugExNetInitSerializer(DrugExNetSerializer):
    molset = serializers.PrimaryKeyRelatedField(many=False, queryset=MolSet.objects.all(), required=False)
    trainingStrategy = DrugExTrainingStrategyInitSerializer(many=False)
    validationStrategy = DrugExValidationStrategyInitSerializer(many=False, required=False)
    parent = serializers.PrimaryKeyRelatedField(many=False, queryset=models.DrugExNet.objects.all(), required=False, allow_null=True)

    class Meta:
        model = models.DrugExNet
        fields = DrugExNetSerializer.Meta.fields
        read_only_fields = DrugExNetSerializer.Meta.read_only_fields

    def create(self, validated_data, **kwargs):
        instance = super().create(
            validated_data,
            molset=validated_data['molset'] if 'molset' in validated_data else None,
            **kwargs
        )
        if "parent" in validated_data and validated_data['parent']:
            instance.parent = validated_data['parent']
            instance.save()

        strat_data = validated_data['trainingStrategy']
        trainingStrategy = models.DrugExNetTraining.objects.create(
            modelInstance=instance,
            algorithm=strat_data['algorithm'],
            mode=strat_data['mode'],
            modelClass=strat_data['modelClass'],
            inputType=strat_data['inputType']
        )

        self.saveParameters(trainingStrategy, strat_data)

        if 'validationStrategy' in validated_data:
            strat_data = validated_data['validationStrategy']
            validationStrategy = models.DrugExNetValidation.objects.create(
                modelInstance = instance,
                validSetSize=strat_data['validSetSize']
            )
            validationStrategy.metrics.set(
                ModelPerformanceMetric.objects
                    .filter(validModes__name=Algorithm.GENERATOR)
                    .filter(Q(validAlgorithms__pk=validated_data['trainingStrategy']['algorithm'].id) | Q(validAlgorithms=None))
                    .distinct()
                    .all()
            )
            validationStrategy.save()

        # create the DrugEx generator with this agent instance
        if instance.molset:
            models.DrugEx.objects.create(
                agent=instance,
                name=instance.name,
                description=instance.description,
                project=instance.project
            )

        return instance


class DrugExAgentTrainingStrategySerializer(TrainingStrategySerializer):

    class Meta:
        model = models.DrugExAgentTraining
        fields = TrainingStrategySerializer.Meta.fields + ('explorer',)


class DrugExAgentTrainingStrategyInitSerializer(TrainingStrategyInitSerializer):

    class Meta:
        model = models.DrugExAgentTraining
        fields = TrainingStrategyInitSerializer.Meta.fields + ('explorer',)


class DrugExAgentValidationStrategySerializer(ValidationStrategySerializer):

    class Meta:
        model = models.DrugExAgentValidation
        fields = ValidationStrategySerializer.Meta.fields + ('validSetSize',)


class DrugExAgentValidationStrategyInitSerializer(ValidationStrategyInitSerializer):
    metrics = serializers.PrimaryKeyRelatedField(many=True, queryset=ModelPerformanceMetric.objects.all(), required=False)

    class Meta:
        model = models.DrugExAgentValidation
        fields = ValidationStrategyInitSerializer.Meta.fields  + ('validSetSize',)

class ScoringFunctionSerializer(serializers.HyperlinkedModelSerializer):
    project = serializers.PrimaryKeyRelatedField(many=False, queryset=Project.objects.all())

    class Meta:
        model = models.ScoringMethod
        fields = ('id', 'name', 'description', 'created', 'updated', 'project')
        read_only_fields = ('id', 'created', 'updated', )

class QSARScorerSerializer(ScoringFunctionSerializer):
    model = serializers.PrimaryKeyRelatedField(many=False, queryset=QSARModel.objects.all())

    class Meta:
        model = models.GenUIModelScorer
        fields = ScoringFunctionSerializer.Meta.fields + ('model',)
        read_only_fields = ('id', 'created', 'updated', )

class PropertyScorerSerializer(ScoringFunctionSerializer):

    class Meta:
        model = models.PropertyScorer
        fields = ScoringFunctionSerializer.Meta.fields + ('prop',)
        read_only_fields = ('id', 'created', 'updated', )

class DrugExScorerSerializer(serializers.HyperlinkedModelSerializer):
    project = serializers.PrimaryKeyRelatedField(many=False, queryset=Project.objects.all())
    environment = serializers.PrimaryKeyRelatedField(many=False, queryset=models.DrugExEnvironment.objects.all())
    modifier = serializers.PrimaryKeyRelatedField(many=False, queryset=models.ScoreModifier.objects.all())
    method = serializers.PrimaryKeyRelatedField(many=False, queryset=models.ScoringMethod.objects.all())

    class Meta:
        model = models.DrugExScorer
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'environment', 'modifier', 'method', 'threshold')
        read_only_fields = ('id', 'created', 'updated', )

class ModifierSerializer(serializers.HyperlinkedModelSerializer):
    project = serializers.PrimaryKeyRelatedField(many=False, queryset=Project.objects.all())

    class Meta:
        model = models.ScoreModifier
        fields = ('id', 'name', 'description', 'created', 'updated', 'project')
        read_only_fields = ('id', 'created', 'updated', )

class ClippedSerializer(ModifierSerializer):

    class Meta:
        model = models.ClippedScore
        fields = ModifierSerializer.Meta.fields + ('upper', 'lower', 'high', 'low', 'smooth')
        read_only_fields = ModifierSerializer.Meta.read_only_fields

class SmoothHumpSerializer(ModifierSerializer):

    class Meta:
        model = models.SmoothHump
        fields = ModifierSerializer.Meta.fields + ('upper', 'lower', 'sigma')
        read_only_fields = ModifierSerializer.Meta.read_only_fields

class DrugExEnvironmentSerializer(serializers.HyperlinkedModelSerializer):
    project = serializers.PrimaryKeyRelatedField(many=False, queryset=Project.objects.all())
    scorers = DrugExScorerSerializer(many=True, read_only=True)

    class Meta:
        model = models.DrugExEnvironment
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'rewardScheme', 'scorers')
        read_only_fields = ('id', 'created', 'updated', 'scorers')

class DrugExEnvironmentCalculationSerializer(serializers.Serializer):
    molsets = serializers.ListField(child=serializers.IntegerField(min_value=1, allow_null=False))
    useModifiers = serializers.BooleanField(allow_null=False, default=True)
    task = serializers.CharField(min_length=1, allow_blank=False, allow_null=True, read_only=True)

    class Meta:
        fields = ('molsets', 'useModifiers', 'task')
        read_only_fields = ('task',)

class DrugExAgentSerializer(ModelSerializer):
    trainingStrategy = DrugExAgentTrainingStrategySerializer(many=False)
    validationStrategy = DrugExAgentValidationStrategySerializer(many=False)
    environment = DrugExEnvironmentSerializer(many=False)
    explorationNet = DrugExNetSerializer(many=False)
    exploitationNet = DrugExNetSerializer(many=False)

    class Meta:
        model = models.DrugExAgent
        fields = [x for x in ModelSerializer.Meta.fields if x not in ('performance',)]  + ["environment", "explorationNet", "exploitationNet"]
        read_only_fields = [x for x in ModelSerializer.Meta.read_only_fields if x not in ('performance',)]


class DrugExAgentInitSerializer(DrugExAgentSerializer):
    trainingStrategy = DrugExAgentTrainingStrategyInitSerializer(many=False)
    validationStrategy = DrugExAgentValidationStrategyInitSerializer(many=False, required=False)
    environment = serializers.PrimaryKeyRelatedField(many=False, queryset=models.DrugExEnvironment.objects.all())
    explorationNet = serializers.PrimaryKeyRelatedField(many=False, queryset=models.DrugExNet.objects.all())
    exploitationNet = serializers.PrimaryKeyRelatedField(many=False, queryset=models.DrugExNet.objects.all())

    class Meta:
        model = models.DrugExAgent
        fields = DrugExAgentSerializer.Meta.fields
        read_only_fields = DrugExAgentSerializer.Meta.read_only_fields

    def create(self, validated_data, **kwargs):
        instance = super().create(
            validated_data,
            environment=validated_data['environment'],
            exploitationNet=validated_data['exploitationNet'],
            explorationNet=validated_data['explorationNet'],
            **kwargs
        )

        strat_data = validated_data['trainingStrategy']
        trainingStrategy = models.DrugExAgentTraining.objects.create(
            modelInstance=instance,
            algorithm=strat_data['algorithm'],
            mode=strat_data['mode'],
            explorer=strat_data['explorer']
        )
        self.saveParameters(trainingStrategy, strat_data)

        if 'validationStrategy' in validated_data:
            strat_data = validated_data['validationStrategy']
            validationStrategy = models.DrugExAgentValidation.objects.create(
                modelInstance = instance
            )
            validationStrategy.validSetSize = strat_data['validSetSize']
            strat_data = {
                'metrics' : ModelPerformanceMetric.objects
                                .filter(validModes__name=Algorithm.GENERATOR)
                                .filter(Q(validAlgorithms__pk=validated_data['trainingStrategy']['algorithm'].id) | Q(validAlgorithms=None))
                                .distinct()
                                .all()
            }
            validationStrategy.metrics.set(strat_data['metrics'])
            validationStrategy.save()

        # create the DrugEx generator with this agent instance
        models.DrugEx.objects.create(
            agent=instance,
            name=instance.name,
            description=instance.description,
            project=instance.project
        )

        return instance


class DrugExGeneratorSerializer(GeneratorSerializer):
    agent = serializers.PrimaryKeyRelatedField(many=False, queryset=Model.objects.filter(trainingStrategies__algorithm__validModes__name__in=('generator',)), required=True)
    molset = serializers.PrimaryKeyRelatedField(many=False, queryset=MolSet.objects.all(), required=False)
    project = serializers.PrimaryKeyRelatedField(many=False, queryset=Project.objects.all(), required=True)
    compounds = serializers.PrimaryKeyRelatedField(many=True, queryset=MolSet.objects.all(), required=False)

    class Meta:
        model = models.DrugEx
        fields = GeneratorSerializer.Meta.fields + ('agent', 'molset')

class ModifierTestSerializer(serializers.Serializer):
    inputs = serializers.ListField(child=serializers.FloatField())
    params = serializers.DictField()
    results = serializers.ListField(child=serializers.FloatField(), required=False, read_only=True)

    class Meta:
        fields = ('inputs', 'params', 'results')
        read_only_fields = ('results',)