from django.db import models

from compounds.models import MolSet, Molecule
from modelling.models import Model, TrainingStrategy
from qsar.models import DescriptorGroup


class Map(Model):
    molsets = models.ManyToManyField(MolSet, related_name="maps")

class Point(models.Model):
    map = models.ForeignKey(Map, on_delete=models.CASCADE, null=False)
    molecule = models.ForeignKey(Molecule, on_delete=models.CASCADE, null=False)
    x = models.FloatField(blank=False, null=False)
    y = models.FloatField(blank=False, null=False)

    class Meta:
        unique_together = ('map', 'molecule',)

class MappingStrategy(TrainingStrategy):
    descriptors = models.ManyToManyField(DescriptorGroup)

