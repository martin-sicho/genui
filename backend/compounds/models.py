from django.db import models
from polymorphic.models import PolymorphicModel
from projects.models import DataSet

# Create your models here.

class MolSet(DataSet):
    pass

class ActivitySet(DataSet):
    pass

class Molecule(PolymorphicModel):

    canonicalSMILES = models.CharField(max_length=65536)
    inchiKey = models.CharField(max_length=65536, unique=True)
    providers = models.ManyToManyField(MolSet, blank=False)

class ChEMBLMolecule(Molecule):
    chemblID = models.CharField(max_length=32, unique=True)

class ChEMBLAssay(models.Model):
    assayID = models.CharField(max_length=32, unique=True)

class ChEMBLTarget(models.Model):
    targetID = models.CharField(max_length=32, unique=True)

class ChEMBLCompounds(MolSet):
    assays = models.ManyToManyField(ChEMBLAssay, blank=False)
    targets = models.ManyToManyField(ChEMBLTarget, blank=False)

class ChEMBLActivities(ActivitySet):
    assays = models.ManyToManyField(ChEMBLAssay, blank=False)
    targets = models.ManyToManyField(ChEMBLTarget, blank=False)

class ActivityUnit(models.Model):
    value = models.CharField(blank=False, max_length=8, unique=True)

class Activity(models.Model):

    value = models.FloatField(blank=False)
    units = models.ForeignKey(ActivityUnit, on_delete=models.CASCADE, blank=False)
    source = models.ForeignKey(ActivitySet, on_delete=models.CASCADE, blank=False)
    molecule = models.ForeignKey(Molecule, on_delete=models.CASCADE, blank=False)

