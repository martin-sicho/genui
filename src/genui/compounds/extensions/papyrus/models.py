from django.db import models
from django.db.models import Count

from genui.compounds.models import Molecule, MolSet, Activity, ActivitySet

class PapyrusAssay(models.Model):
    assayID = models.CharField(max_length=2800, unique=True, blank=False) # Longest string found in the dataset is 2499 characters

class PapyrusTarget(models.Model):
    targetID = models.CharField(max_length=128, unique=True, blank=False)

class PapyrusMolecule(Molecule):
    CID = models.CharField(max_length=2800, unique=False, blank=True, null=False)

class PapyrusCompounds(MolSet):
    targets = models.ManyToManyField(PapyrusTarget, blank=False)

class PapyrusActivity(Activity):
    relation = models.CharField(blank=False, max_length=128)
    assay = models.ForeignKey(PapyrusAssay, on_delete=models.CASCADE, null=False, blank=False)
    target = models.ForeignKey(PapyrusTarget, on_delete=models.CASCADE, null=False, blank=False)
    comment = models.CharField(blank=True, max_length=128, null=True)

class PapyrusActivities(ActivitySet):

    class PapyrusActivitySetSummary(ActivitySet.ActivitySetSummary):

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
        return self.PapyrusActivitySetSummary(self)

    def cleanForModelling(self, activity_type):
        activities = []
        mols = []
        units = None
        for activity in PapyrusActivity.objects.filter(source=self, type=activity_type, relation="="):
            units = activity.units
            mols.append(activity.molecule)
            activities.append(activity.value)

        return mols, activities, units
