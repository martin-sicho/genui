from django.db import models

from compounds.models import MolSet, Molecule
from modelling.models import Model, TrainingStrategy
from qsar.models import DescriptorGroup


class Map(Model):
    molsets = models.ManyToManyField(MolSet, related_name="maps")

class Point(models.Model):
    molecule = models.ForeignKey(Molecule, on_delete=models.CASCADE, null=False)
    x = models.FloatField(blank=False, null=False)
    y = models.FloatField(blank=False, null=False)

class MappingStrategy(TrainingStrategy):
    descriptors = models.ManyToManyField(DescriptorGroup)

