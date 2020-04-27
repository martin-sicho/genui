from django.db.models import Count
from django_rdkit import models
from django_celery_results.models import TaskResult
from djcelery_model.models import TaskMixin
from polymorphic.models import PolymorphicModel
from rdkit import Chem

from genui.commons.models import TaskShortcutsMixIn, PolymorphicTaskManager, OverwriteStorage
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

    molecules = models.ForeignKey(MolSet, blank=False, null=True, on_delete=models.CASCADE, related_name="activities") # FIXME: it probably makes more sense to make this field non-nullable

    # TODO: arguments should be added to this method that allow specification of the activity type and units
    def cleanForModelling(self, activity_type : "ActivityTypes") -> tuple:
        """
        All subclasses should override this method to implement a procedure that returns
        molecules as Molecule instances and their activities ready for modelling.

        :return: Tuple of list objects (same length) -> Molecule instances and their associated activity values for modelling
        """
        raise NotImplementedError("This should be overridden in children, which seems not to be the case.")

    def getSummary(self):
        return self.ActivitySetSummary(self)

class Molecule(PolymorphicModel):
    canonicalSMILES = models.CharField(max_length=65536, unique=True, blank=False)
    inchiKey = models.CharField(max_length=65536, unique=True, blank=False)
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

    extension = models.CharField(choices=FILE_TYPES, null=False, blank=False, default=SVG, max_length=128)

class MoleculePic(models.Model):
    format = models.ForeignKey(PictureFormat, on_delete=models.CASCADE)
    molecule = models.ForeignKey(Molecule, on_delete=models.CASCADE, related_name='pics')
    image = models.ImageField(upload_to='compounds/pics/', storage=OverwriteStorage())

class ChEMBLAssay(models.Model):
    assayID = models.CharField(max_length=128, unique=True, blank=False)

class ChEMBLTarget(models.Model):
    targetID = models.CharField(max_length=128, unique=True, blank=False)

class ChEMBLMolecule(Molecule):
    chemblID = models.CharField(max_length=128, unique=True, blank=False, null=False)

class ChEMBLCompounds(MolSet):
    targets = models.ManyToManyField(ChEMBLTarget, blank=False)

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
    source = models.ForeignKey(ActivitySet, on_delete=models.CASCADE, blank=False, related_name='activities')
    molecule = models.ForeignKey(Molecule, on_delete=models.CASCADE, blank=False, related_name='activities')
    parent = models.ForeignKey("self", null=True, on_delete=models.CASCADE, related_name="children")

    def __str__(self):
        return '%s object (%s=%f)' % (self.__class__.__name__, self.type.value, self.value)

class ChEMBLActivity(Activity):
    relation = models.CharField(blank=False, max_length=128)
    assay = models.ForeignKey(ChEMBLAssay, on_delete=models.CASCADE, null=False, blank=False)
    target = models.ForeignKey(ChEMBLTarget, on_delete=models.CASCADE, null=False, blank=False)
    comment = models.CharField(blank=True, max_length=128, null=True)

class ChEMBLActivities(ActivitySet):

    class ChEMBLActivitySetSummary(ActivitySet.ActivitySetSummary):

        def getTypeInfo(self):
            # return self.activities.filter(chemblactivity__relation="=").values('type').annotate(
            #     occurences=Count('id'),
            #     molecules=Count('molecule', distinct=True),
            # ).order_by('-molecules')
            return self.activities.values('type').annotate(
                    occurences=Count('id'),
                    molecules=Count('molecule', distinct=True),
                ).order_by('-molecules')

    def getSummary(self):
        return self.ChEMBLActivitySetSummary(self)

    def cleanForModelling(self, activity_type):
        # FIXME: we could do much more here -> remove weird molecules and stuff (see some good QSAR preprocessing workflow)
        activities = []
        mols = []
        units = None
        for activity in ChEMBLActivity.objects.filter(source=self, type=activity_type, relation="="):
            units = activity.units # FIXME: make sure that all activities are in the same units
            mols.append(activity.molecule)
            activities.append(activity.value)

        return mols, activities, units
