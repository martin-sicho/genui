from django.db.models import Count
from django_rdkit import models
from djcelery_model.models import TaskMixin
from polymorphic.models import PolymorphicModel
from rdkit import Chem
from rdkit.Chem import AllChem

from genui.utils.models import OverwriteStorage, NON_POLYMORPHIC_CASCADE
from genui.utils.extensions.tasks.models import TaskShortcutsMixIn, PolymorphicTaskManager
from . import helpers
from genui.projects.models import DataSet

class MolSet(TaskShortcutsMixIn, TaskMixin, DataSet):
    objects = PolymorphicTaskManager()

    def __str__(self):
        return '%s object (%s)' % (self.__class__.__name__, self.name)

class ActivitySet(TaskShortcutsMixIn, TaskMixin, DataSet):
    objects = PolymorphicTaskManager()

    class ActivityTypeSummary:

        def __init__(self, activityType : "ActivityTypes", compounds : "Molecule", occurences: "Activity"):
            self.type = activityType
            self.moleculesTotal = compounds
            self.activitiesTotal = occurences

    class ActivitySetSummary:

        def generateTypeSummaries(self):
            typeInfo = self.getTypeInfo()
            ret = []
            for info in typeInfo:
                actype = ActivityTypes.objects.get(pk=info['type'])
                ret.append(ActivitySet.ActivityTypeSummary(
                    actype,
                    info['molecules'],
                    info['occurences']
                ))
            return ret

        def getTypeInfo(self):
            return self.activities.values('type').annotate(
                occurences=Count('id'),
                molecules=Count('molecule', distinct=True)
            ).order_by('-molecules')

        def __init__(self, activitySet : "ActivitySet"):
            self._activitySet = activitySet
            self.activities = Activity.objects.filter(source=self._activitySet)
            self.molecules = Molecule.objects.filter(activities__source=self._activitySet).distinct()
            self.typeSummaries = self.generateTypeSummaries()
            self.activitiesTotal = self.activities.count()
            self.moleculesTotal = self.molecules.count()

    molecules = models.ForeignKey(MolSet, blank=False, null=False, on_delete=models.CASCADE, related_name="activities")

    def cleanForModelling(self, activity_type : "ActivityTypes") -> tuple:
        """
        All subclasses should override this method to implement a procedure that returns
        molecules as Molecule instances and their activities ready for models.

        :return: Tuple of list objects (same length) -> Molecule instances and their associated activity values for models
        """
        raise NotImplementedError("This should be overridden in children, which seems not to be the case.")

    def getSummary(self):
        return self.ActivitySetSummary(self)

class ChemicalEntity(models.Model):
    canonicalSMILES = models.CharField(max_length=65536, unique=True, blank=False)
    inchi = models.CharField(max_length=65536, unique=True, blank=False)
    inchiKey = models.CharField(max_length=65536, unique=True, blank=False)
    # from django-rdkit
    rdMol = models.MolField()
    morganFP = models.BfpField(null=True)

    class Meta:
        unique_together = ('canonicalSMILES', 'inchiKey')

    def __str__(self):
        return '%s object <%s>' % (self.__class__.__name__, self.inchiKey)

    @property
    def fingerprint(self):
        if not self.morganFP:
            self.morganFP = AllChem.GetMorganFingerprintAsBitVect(self.rdMol, radius=2, nBits=512)
        return self.morganFP

class Molecule(PolymorphicModel):
    entity = models.ForeignKey(ChemicalEntity, on_delete=models.CASCADE, null=False, related_name='molecules')
    providers = models.ManyToManyField(MolSet, blank=False, related_name='molecules')

    @classmethod
    def create(cls, canonicalSMILES, inchiKey, rdMol, *args, **kwargs):
        ret = cls.objects.create(
            entity=ChemicalEntity.objects.get_or_create(
                canonicalSMILES=canonicalSMILES,
                inchiKey=inchiKey,
                rdMol=rdMol
            )[0],
            *args,
            **kwargs
        )
        return ret

    def __str__(self):
        return '%s object <%s>' % (self.__class__.__name__, self.smiles)

    @property
    def canonicalSMILES(self):
        return self.entity.canonicalSMILES

    @property
    def inchi(self):
        return self.entity.inchi

    @property
    def inchiKey(self):
        return self.entity.inchiKey

    @property
    def rdMol(self):
        return self.entity.rdMol

    @property
    def fingerprint(self):
        return self.entity.fingerprint

    @property
    def smiles(self):
        """
        A shorthand to get a nice human readable SMILES string directly from representation.

        """

        return Chem.MolToSmiles(self.rdMol, isomericSmiles=True, canonical=True)

    def getPic(self, format):
        format = PictureFormat.objects.get_or_create(extension=f'.{format}')[0]
        qs = self.pics.filter(format=format)
        pic = qs.all()[0] if qs.exists() else helpers.createPic(self, format)
        return pic

    @property
    def mainPic(self):
        return self.getPic('svg')

class PictureFormat(models.Model):
    PNG = ".png"
    SVG = ".svg"
    FILE_TYPES = [
       (PNG, 'PNG'),
       (SVG, 'SVG'),
    ]

    extension = models.CharField(choices=FILE_TYPES, null=False, blank=False, default=SVG, max_length=128, unique=True)

class MoleculePic(models.Model):
    format = models.ForeignKey(PictureFormat, on_delete=models.CASCADE)
    molecule = models.ForeignKey(Molecule, on_delete=models.CASCADE, related_name='pics')
    image = models.ImageField(upload_to='compounds/pics/', storage=OverwriteStorage())

class ActivityUnits(models.Model):
    value = models.CharField(blank=False, max_length=128, unique=True)

    def __str__(self):
        return '%s object <%s>' % (self.__class__.__name__, self.value)

class ActivityTypes(models.Model):
    value = models.CharField(blank=False, max_length=128, unique=True)

    def __str__(self):
        return '%s object <%s>' % (self.__class__.__name__, self.value)

class Activity(PolymorphicModel):
    value = models.FloatField(blank=False)
    type = models.ForeignKey(ActivityTypes, on_delete=models.CASCADE, null=False)
    units = models.ForeignKey(ActivityUnits, on_delete=models.CASCADE, null=True)
    source = models.ForeignKey(ActivitySet, on_delete=NON_POLYMORPHIC_CASCADE, blank=False, related_name='activities')
    molecule = models.ForeignKey(Molecule, on_delete=models.CASCADE, blank=False, related_name='activities')
    parent = models.ForeignKey("self", null=True, on_delete=models.CASCADE, related_name="children")

    def __str__(self):
        return '%s object (%s=%f)' % (self.__class__.__name__, self.type.value, self.value)

class MolSetFile(PolymorphicModel):
    molset = models.ForeignKey(MolSet, on_delete=models.CASCADE, null=False, related_name='files')
    file = models.FileField(null=False, upload_to='compounds/sets/files/', storage=OverwriteStorage())

    @staticmethod
    def create(molset, filename, file):
        file.name = filename
        return MolSetFile.objects.create(molset=molset, file=file)
