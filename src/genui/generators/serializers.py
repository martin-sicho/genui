"""
serializers

Created by: Martin Sicho
On: 27-01-20, 17:00
"""
from django.db.models import Q
from rest_framework import serializers

from genui.utils.serializers import GenericModelSerializerMixIn
from genui.compounds.models import MolSet
from genui.compounds.serializers import MolSetSerializer, GenericMolSetSerializer
from genui.modelling.core.bases import Algorithm
from genui.modelling.models import ModelPerformanceMetric
from genui.modelling.serializers import ModelSerializer, ValidationStrategySerializer, TrainingStrategySerializer, \
    TrainingStrategyInitSerializer, ValidationStrategyInitSerializer
from genui.projects.serializers import ProjectSerializer
from genui.qsar.models import QSARModel
from genui.qsar.serializers import QSARModelSerializer
from . import models

class GeneratorSerializer(GenericModelSerializerMixIn, serializers.HyperlinkedModelSerializer):
    className = GenericModelSerializerMixIn.className
    extraArgs = GenericModelSerializerMixIn.extraArgs
    project = ProjectSerializer(many=False)
    compounds = MolSetSerializer(many=True)

    class Meta:
        model = models.Generator
        fields = ('id', 'name', 'description', 'project', 'compounds', 'className', 'extraArgs')

class GeneratedSetSerializer(GenericMolSetSerializer):
    source = GeneratorSerializer(many=False)

    class Meta:
        model = models.GeneratedMolSet
        fields = list(GenericMolSetSerializer.Meta.fields) + ['source']
        read_only_fields = list(GenericMolSetSerializer.Meta.read_only_fields)

class GeneratedSetInitSerializer(GeneratedSetSerializer):
    source = serializers.PrimaryKeyRelatedField(many=False, queryset=models.Generator.objects.all())
    nSamples = serializers.IntegerField(min_value=1)
    taskID = serializers.UUIDField(required=False, read_only=True)

    class Meta:
        model = models.GeneratedMolSet
        fields = GeneratedSetSerializer.Meta.fields + ["nSamples", "taskID"]
        read_only_fields = GeneratedSetSerializer.Meta.read_only_fields + ["taskID"]

    def create(self, validated_data):
        molset_class = models.GeneratedMolSet
        if "className" in validated_data:
            # TODO:  make this more general (look in different modules as well)
            molset_class = getattr(models, validated_data["className"])
        extraArgs = validated_data["extraArgs"] if "extraArgs" in validated_data else dict()
        instance = molset_class.objects.create(
            name=validated_data["name"],
            description=validated_data["description"] if "description" in validated_data else "",
            source=validated_data["source"],
            project=validated_data["project"],
            **extraArgs

        )
        instance.nSamples = validated_data["nSamples"]

        return instance

class DrugExValidationStrategySerializer(ValidationStrategySerializer):

    class Meta:
        model = models.DrugExValidationStrategy
        fields = ValidationStrategySerializer.Meta.fields + ("validSetSize",)

class DrugExValidationStrategyInitSerializer(DrugExValidationStrategySerializer):
    metrics = serializers.PrimaryKeyRelatedField(many=True, queryset=ModelPerformanceMetric.objects.all(), required=False)

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
        trainingStrategy = models.DrugExNetTrainingStrategy.objects.create(
            modelInstance=instance,
            algorithm=strat_data['algorithm'],
            mode=strat_data['mode'],
        )

        self.saveParameters(trainingStrategy, strat_data)

        if 'validationStrategy' in validated_data:
            strat_data = validated_data['validationStrategy']
            validationStrategy = models.DrugExValidationStrategy.objects.create(
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
        models.DrugEx.objects.create(
            agent=instance,
            name=instance.name,
            description=instance.description,
            project=instance.project
        )

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

class DrugExAgentValidationStrategyInitSerializer(ValidationStrategyInitSerializer):

    class Meta:
        model = models.DrugExAgentValidationStrategy
        fields = ValidationStrategyInitSerializer.Meta.fields

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

    def create(self, validated_data, **kwargs):
        instance = super().create(
            validated_data,
            environment=validated_data['environment'],
            exploitationNet=validated_data['exploitationNet'],
            explorationNet=validated_data['explorationNet'],
            **kwargs
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
    agent = DrugExAgentSerializer(many=False)

    class Meta:
        model = models.DrugEx
        fields = GeneratorSerializer.Meta.fields + ('agent',)