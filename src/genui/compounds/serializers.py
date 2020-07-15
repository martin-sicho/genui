"""
serializers

Created by: Martin Sicho
On: 18-12-19, 10:27
"""
from rest_framework import serializers

from genui.utils.serializers import GenericModelSerializerMixIn
from genui.projects.models import Project
from .models import MolSet, Molecule, MoleculePic, PictureFormat, \
    ActivitySet, Activity, ActivityUnits, ActivityTypes, MolSetFile


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
    # pics = MoleculePicSerializer(many=True, required=False)
    mainPic = MoleculePicSerializer(many=False, required=True)
    properties = serializers.SerializerMethodField(required=False)

    class Meta:
        model = Molecule
        fields = ('id', 'smiles', 'inchiKey', 'providers', 'mainPic',  'properties', 'className', 'extraArgs')

    def get_properties(self, obj):
        props = [x for x in dir(obj) if x.startswith('rdkit_prop_')]
        return {prop.split('_')[-1] : getattr(obj, prop) for prop in props}

class MolSetFileSerializer(serializers.HyperlinkedModelSerializer):
    molset = serializers.PrimaryKeyRelatedField(queryset=MolSetFile.objects.all(), required=True)

    class Meta:
        model = MolSetFile
        fields = ('id', 'molset', 'file')
        read_only_fields = ('id',)

    def create(self, validated_data):
        return MolSetFile.create(
            molset=validated_data['molset'],
            filename=validated_data['file'].name,
            file=validated_data['file'],
        )

class MolSetSerializer(serializers.HyperlinkedModelSerializer):
    project = serializers.PrimaryKeyRelatedField(many=False, queryset=Project.objects.all())
    activities = serializers.PrimaryKeyRelatedField(many=True, read_only=True)
    files = MolSetFileSerializer(many=True, required=False, allow_null=False, read_only=True)

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
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'activities', 'files')
        read_only_fields = ('created', 'updated', 'activities', 'files')

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

class ActivityUnitsSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = ActivityUnits
        fields = ('id', 'value',)

class ActivityTypeSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = ActivityTypes
        fields = ('id', 'value',)

class ActivityTypeSummary(serializers.Serializer):
    type = ActivityTypeSerializer(many=False)
    moleculesTotal = serializers.IntegerField(min_value=0, required=True)
    activitiesTotal = serializers.IntegerField(min_value=0, required=True)

class ActivitySetSummarySerializer(serializers.Serializer):
    moleculesTotal = serializers.IntegerField(min_value=0, required=True)
    activitiesTotal = serializers.IntegerField(min_value=0, required=True)
    typeSummaries = ActivityTypeSummary(many=True)

class ActivitySerializer(GenericModelSerializerMixIn, serializers.HyperlinkedModelSerializer):
    className = GenericModelSerializerMixIn.className
    extraArgs = GenericModelSerializerMixIn.extraArgs

    units = ActivityUnitsSerializer(many=False)
    type = ActivityTypeSerializer(many=False)
    source = serializers.PrimaryKeyRelatedField(many=False, queryset=ActivitySet.objects.all())
    molecule = serializers.PrimaryKeyRelatedField(many=False, queryset=Molecule.objects.all())
    parent = serializers.PrimaryKeyRelatedField(many=False, queryset=Activity.objects.all())

    class Meta:
        model = Activity
        fields = ('id', 'value', 'type', 'units', 'source', 'molecule', 'parent', 'className', 'extraArgs')