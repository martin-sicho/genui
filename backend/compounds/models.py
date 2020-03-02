from django_rdkit import models
from django_celery_results.models import TaskResult
from djcelery_model.models import TaskMixin
from polymorphic.models import PolymorphicModel
from rdkit import Chem

from commons.models import TaskShortcutsMixIn, PolymorphicTaskManager
from projects.models import DataSet

class MolSet(TaskShortcutsMixIn, TaskMixin, DataSet):
    objects = PolymorphicTaskManager()

    def __str__(self):
        return '%s object (%s)' % (self.__class__.__name__, self.name)

class ActivitySet(TaskShortcutsMixIn, TaskMixin, DataSet):
    objects = PolymorphicTaskManager()

    molecules = models.ForeignKey(MolSet, blank=False, null=True, on_delete=models.CASCADE, related_name="activities") # FIXME: it probably makes more sense to make this field non-nullable

class Molecule(PolymorphicModel):
    canonicalSMILES = models.CharField(max_length=65536, unique=True)
    inchiKey = models.CharField(max_length=65536, unique=True)
    providers = models.ManyToManyField(MolSet, blank=False, related_name='molecules')

    # from django-rdkit
    molObject = models.MolField()
    morganFP2 = models.BfpField(null=True)

    def __str__(self):
        return '%s object <%s>' % (self.__class__.__name__, self.smiles)

    @property
    def smiles(self):
        """
        A shorthand to get a nice human readable SMILES string directly from representation.

        """

        return Chem.MolToSmiles(self.molObject, isomericSmiles=True, canonical=True)

class PictureFormat(models.Model):
    PNG = ".png"
    SVG = ".svg"
    FILE_TYPES = [
       (PNG, 'PNG'),
       (SVG, 'SVG'),
    ]

    extension = models.CharField(choices=FILE_TYPES, null=False, blank=False, default=SVG, max_length=16)

class MoleculePic(models.Model):
    format = models.ForeignKey(PictureFormat, on_delete=models.CASCADE)
    molecule = models.ForeignKey(Molecule, on_delete=models.CASCADE, related_name='pics')
    image = models.ImageField(upload_to='compounds/pics/')

class ChEMBLAssay(models.Model):
    assayID = models.CharField(max_length=32, unique=True, blank=False)

class ChEMBLTarget(models.Model):
    targetID = models.CharField(max_length=32, unique=True, blank=False)

class ChEMBLMolecule(Molecule):
    chemblID = models.CharField(max_length=32, unique=True, blank=False, null=False)

class ChEMBLCompounds(MolSet):
    targets = models.ManyToManyField(ChEMBLTarget, blank=False)

class ActivityUnits(models.Model):
    value = models.CharField(blank=False, max_length=8, unique=True)

class Activity(models.Model):
    value = models.FloatField(blank=False)
    units = models.ForeignKey(ActivityUnits, on_delete=models.CASCADE, null=True)
    source = models.ForeignKey(ActivitySet, on_delete=models.CASCADE, blank=False, related_name='activities')
    molecule = models.ForeignKey(Molecule, on_delete=models.CASCADE, blank=False, related_name="activities")

class ChEMBLActivity(Activity):
    type = models.CharField(blank=False, max_length=128)
    relation = models.CharField(blank=False, max_length=128)
    assay = models.ForeignKey(ChEMBLAssay, on_delete=models.CASCADE, null=False, blank=False)
    target = models.ForeignKey(ChEMBLTarget, on_delete=models.CASCADE, null=False, blank=False)
    comment = models.CharField(blank=True, max_length=128, null=True)

class ChEMBLActivities(ActivitySet):

    def cleanForModelling(self):
        activities = []
        mols = []
        for activity in ChEMBLActivity.objects.filter(source=self):
            if activity.type == "PCHEMBL_VALUE" and activity.relation == "=":
                mols.append(activity.molecule)
                activities.append(activity.value)

        return mols, activities
