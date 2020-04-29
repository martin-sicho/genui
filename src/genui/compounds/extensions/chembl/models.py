from django.db import models
from django.db.models import Count

from genui.compounds.models import Molecule, MolSet, Activity, ActivitySet


class ChEMBLAssay(models.Model):
    assayID = models.CharField(max_length=128, unique=True, blank=False)

class ChEMBLTarget(models.Model):
    targetID = models.CharField(max_length=128, unique=True, blank=False)

class ChEMBLMolecule(Molecule):
    chemblID = models.CharField(max_length=128, unique=True, blank=False, null=False)

class ChEMBLCompounds(MolSet):
    targets = models.ManyToManyField(ChEMBLTarget, blank=False)

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
