"""
serializers

Created by: Martin Sicho
On: 18-12-19, 10:27
"""
from rest_framework import serializers

from commons.serializers import GenericModelSerializerMixIn
from projects.models import Project
from .models import MolSet, Molecule, ChEMBLCompounds, ChEMBLTarget, ChEMBLAssay, MoleculePic, PictureFormat, \
    ActivitySet, Activity, ActivityUnits, ActivityTypes


class PictureFormatSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = PictureFormat
        fields = ('id', 'extension',)

class MoleculePicSerializer(serializers.HyperlinkedModelSerializer):
    format = PictureFormatSerializer(many=False)
    molecule = serializers.PrimaryKeyRelatedField(queryset=Molecule.objects.all(), many=False)

    class Meta:
        model = MoleculePic
        fields = ('format', 'image', 'molecule')

class MoleculeSerializer(GenericModelSerializerMixIn, serializers.HyperlinkedModelSerializer):
    className = GenericModelSerializerMixIn.className
    extraArgs = GenericModelSerializerMixIn.extraArgs

    providers = serializers.PrimaryKeyRelatedField(many=True, queryset=MolSet.objects.all())
    pics = MoleculePicSerializer(many=True, required=False)
    mainPic = MoleculePicSerializer(many=False, required=True)

    class Meta:
        model = Molecule
        fields = ('id', 'smiles', 'inchiKey', 'providers', 'mainPic', 'pics',  'className', 'extraArgs')

class ChEMBLAssaySerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = ChEMBLAssay
        fields = ('assayID',)

class ChEMBLTargetSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = ChEMBLTarget
        fields = ('targetID',)

class MolSetSerializer(serializers.HyperlinkedModelSerializer):
    project = serializers.PrimaryKeyRelatedField(many=False, queryset=Project.objects.all())
    activities = serializers.PrimaryKeyRelatedField(many=True, queryset=ActivitySet.objects.all())

    class AutoSchemaMixIn:
        def get_operation(self, path, method):
            ret = super().get_operation(path, method)
            if method in ('POST', 'PUT', 'PATCH'):
                ret['responses']['200']['content']['application/json']['schema']['properties']['task'] = {
                    'type' : 'string'
                }
            return ret

    class Meta:
        model = MolSet
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'activities')
        read_only_fields = ('created', 'updated', 'activities')

class GenericMolSetSerializer(GenericModelSerializerMixIn, MolSetSerializer):
    className = GenericModelSerializerMixIn.className
    extraArgs = GenericModelSerializerMixIn.extraArgs

    class Meta:
        model = MolSet
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'activities', 'className', 'extraArgs')
        read_only_fields = ('created', 'updated', 'extraArgs', 'activities')

class ActivitySetSerializer(GenericModelSerializerMixIn, serializers.HyperlinkedModelSerializer):
    className = GenericModelSerializerMixIn.className
    extraArgs = GenericModelSerializerMixIn.extraArgs

    project = serializers.PrimaryKeyRelatedField(many=False, queryset=Project.objects.all())
    molecules = serializers.PrimaryKeyRelatedField(many=False, queryset=MolSet.objects.all())

    class Meta:
        model = ActivitySet
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'molecules', 'className', 'extraArgs')
        read_only_fields = ('created', 'updated', 'className')

class ActivityUnitSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = ActivityUnits
        fields = ('id', 'value',)

class ActivityTypeSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = ActivityTypes
        fields = ('id', 'value',)

class ActivitySerializer(GenericModelSerializerMixIn, serializers.HyperlinkedModelSerializer):
    className = GenericModelSerializerMixIn.className
    extraArgs = GenericModelSerializerMixIn.extraArgs

    units = ActivityUnitSerializer(many=False)
    type = ActivityTypeSerializer(many=False)
    source = serializers.PrimaryKeyRelatedField(many=False, queryset=ActivitySet.objects.all())
    molecule = serializers.PrimaryKeyRelatedField(many=False, queryset=Molecule.objects.all())

    class Meta:
        model = Activity
        fields = ('id', 'value', 'type', 'units', 'source', 'molecule', 'className', 'extraArgs')

class ChEMBLSetSerializer(MolSetSerializer):
    targets = ChEMBLTargetSerializer(many=True)

    class Meta:
        model = ChEMBLCompounds
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'targets', 'activities')
        read_only_fields = ('created', 'updated')

class ChEMBLSetInitSerializer(ChEMBLSetSerializer):
    maxPerTarget = serializers.IntegerField(min_value=1, required=False)
    taskID = serializers.CharField(required=False, read_only=True)
    targets = serializers.ListField(child=serializers.CharField(), min_length=1, required=True, write_only=True)

    def create(self, validated_data):
        instance = ChEMBLCompounds(
            name=validated_data["name"]
            , description=validated_data["description"]
            , project=validated_data["project"]
        )
        instance.save()
        for target in validated_data['targets']:
            target = ChEMBLTarget.objects.get_or_create(targetID=target)[0]
            instance.targets.add(target)
        instance.save()
        return instance

    class Meta:
        model = ChEMBLCompounds
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'maxPerTarget', 'taskID', 'targets')
        read_only_fields = ('created', 'updated', 'taskID', 'targets')

class ChEMBLSetUpdateSerializer(ChEMBLSetInitSerializer):
    project = serializers.PrimaryKeyRelatedField(many=False, queryset=Project.objects.all(), required=False)
    targets = ChEMBLTargetSerializer(many=True, read_only=True, required=False)
    name = serializers.CharField(required=False, max_length=ChEMBLCompounds._meta.get_field('name').max_length)

    class Meta:
        model = ChEMBLCompounds
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'taskID', 'targets')
        read_only_fields = ('created', 'updated', 'taskID', 'targets')

    def update(self, instance, validated_data):
        for (key, value) in validated_data.items():
            setattr(instance, key, value)
        instance.save()
        return instance
