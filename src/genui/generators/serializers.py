"""
serializers

Created by: Martin Sicho
On: 27-01-20, 17:00
"""
from rest_framework import serializers

from genui.utils.serializers import GenericModelSerializerMixIn
from genui.compounds.serializers import MolSetSerializer
from genui.projects.serializers import ProjectSerializer
from . import models

class GeneratorSerializer(GenericModelSerializerMixIn, serializers.HyperlinkedModelSerializer):
    className = GenericModelSerializerMixIn.className
    extraArgs = GenericModelSerializerMixIn.extraArgs
    project = ProjectSerializer(many=False)
    compounds = MolSetSerializer(many=True)

    class Meta:
        model = models.Generator
        fields = ('id', 'name', 'description', 'project', 'compounds', 'className', 'extraArgs')


