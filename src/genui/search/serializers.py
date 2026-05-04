"""
serializers.py in src/genui/search

"""
from rest_framework import serializers
from genui.compounds.serializers import MoleculeSerializer, ActivitySerializer, MolSetSerializer
from genui.projects.serializers import ProjectSerializer
from genui.compounds.models import Molecule, ChemicalEntity, MolSet, Activity
from genui.projects.models import Project
from rdkit import Chem


class BaseSearchRequestSerializer(serializers.Serializer):
    input = serializers.CharField(required=True, allow_blank=False, trim_whitespace=True, error_messages={"required": "Enter a valid structure"})
    canonical = serializers.CharField(read_only=True)

    ids = serializers.ListField(
        child=serializers.IntegerField(min_value=1, required=True),
        allow_empty=True,
        required=False,
        default=[]
    )

    def parse_input(self, s:str):
        raise NotImplementedError
    
    def canonicalize(self, s:str):
        raise NotImplementedError

    def validate(self, attrs):
        s = attrs.get("input", "")
        mol = self.parse_input(s)
        attrs["canonical"] = self.canonicalize(mol)
        return attrs
    

class SimilaritySearchRequestSerializer(BaseSearchRequestSerializer):

    FP_CHOICES = (
        ("morganFP", "Morgan (ECFP-like)"),
        ("maccsFP", "MACCS Keys"),
    )
    METRIC_CHOICES = (
        ("tanimoto", "Tanimoto"),
        ("dice", "Dice"),
    )

    threshold = serializers.FloatField(min_value=0.0, max_value=1.0, required=False, allow_null=True, default=None)
    fp_type = serializers.ChoiceField(choices=FP_CHOICES, required=False, default="morganFP")
    metric = serializers.ChoiceField(choices=METRIC_CHOICES, required=False, default="tanimoto")
    top_n = serializers.IntegerField(min_value=1, max_value=1000, required=False, default=5)

    def parse_input(self, s: str):
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            raise serializers.ValidationError({"input": "Invalid SMILES."})
        return mol

    def canonicalize(self, mol):
        return Chem.MolToSmiles(mol, canonical=True)

    
    
class SubstructureSearchRequestSerializer(BaseSearchRequestSerializer):

    def parse_input(self, s: str):
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            raise serializers.ValidationError({"input": "Invalid SMILES."})
        return mol

    def canonicalize(self, mol):
        return Chem.MolToSmiles(mol, canonical=True)

    
class SmartsSearchRequestSerializer(BaseSearchRequestSerializer):

    def parse_input(self, s: str):
        mol = Chem.MolFromSmarts(s)
        if mol is None:
            raise serializers.ValidationError({"input": "Invalid SMARTS."})
        return mol
    
    def canonicalize(self, mol):
        return Chem.MolToSmarts(mol)
    

class InchiKeySearchRequestSerializer(serializers.Serializer):
    input = serializers.CharField(required=True, allow_blank=False, trim_whitespace=True, error_messages={"required": "Enter a valid inchiKey"})
    

# class PropertyFilterSerializer(serializers.Serializer):
    
#     PROPERTY_CHOICES = (
#         ("amw", "Molar Weight"),
#         ("hbd", "Hydrogen Bond Donor"),
#         ("hba", "Hydrogen Bond Acceptor"),
#         ("logp", "logP"),
#         ("tpsa","Topological Polar Surface Area")
#     )

#     RELATION_CHOICES = (
#         ("exact","Equals"),
#         ("gt", "Greater Than"),
#         ("gte","Greater Than, Equal"),
#         ("lt","Less Than"),
#         ("lte","Less Than, Equal")
#     )

#     property = serializers.ChoiceField(choices=PROPERTY_CHOICES, required=True)
#     relation = serializers.ChoiceField(choices=RELATION_CHOICES, required=True)
#     value = serializers.FloatField(required=True, allow_null=False)
#     top_n = serializers.IntegerField(min_value=1, max_value=1000, required=False, default=5)

#     ids = serializers.ListField(
#         child=serializers.IntegerField(min_value=1, required=True),
#         allow_empty=True,
#         required=False,
#         default=[]
#     )


class HitSerializer(MoleculeSerializer):
    similarity = serializers.FloatField(read_only=True,required=False)
    project_ids = serializers.ListField(
        child=serializers.IntegerField(read_only=True),
        required=False
    )

    class Meta:
        model = Molecule
        fields = MoleculeSerializer.Meta.fields + ("similarity","project_ids",)


class OccurrenceSerializer(ProjectSerializer):
    providers = serializers.ListField(
         child=serializers.DictField(read_only=True)
     )

    class Meta:
        model = Project
        fields = ("id", "name", "providers")


class BaseSearchResponseSerializer(serializers.Serializer):
    query = serializers.DictField(read_only=True) 

    hits = HitSerializer(many=True, read_only=True)
    total_searched = serializers.IntegerField(read_only=True, required=False)
    total_returned = serializers.IntegerField(read_only=True, required=False)


class SimilaritySearchResponseSerializer(BaseSearchResponseSerializer):
    query = SimilaritySearchRequestSerializer(read_only=True)


class SubstructureSearchResponseSerializer(BaseSearchResponseSerializer):
    query = SubstructureSearchRequestSerializer(read_only=True)


class SmartsSearchResponseSerializer(BaseSearchResponseSerializer):
    query = SmartsSearchRequestSerializer(read_only=True)


class InchiKeySearchResponseSerializer(BaseSearchResponseSerializer):
    query = InchiKeySearchRequestSerializer(read_only=True)


class OccurrenceSearchResponseSerializer(serializers.Serializer):
    query = InchiKeySearchRequestSerializer(read_only=True) 
    occurrence = OccurrenceSerializer(read_only=True, many=True)


# class PropertyFiltersResponseSerializer(BaseSearchResponseSerializer):
#     query = PropertyFilterSerializer(read_only=True)
