from django.db import models

from genui.compounds.extensions.fileimports.model import FileCompounds
from genui.compounds.models import Molecule

class SDFCompounds(FileCompounds):
    activitiesProp = models.CharField(max_length=256, null=False, default='GENUI_ACTIVITIES', blank=False)
    activityTypesProp = models.CharField(max_length=256, null=False, default='GENUI_ACTIVITY_TYPES', blank=False)
    activityUnitsProp = models.CharField(max_length=256, null=True, default='GENUI_ACTIVITY_UNITS', blank=True)
    dataSeparator = models.CharField(max_length=1, null=False, default=',', blank=False)

class SDFMolecule(Molecule):
    name = models.CharField(max_length=1024, null=False, blank=True)

    def __str__(self):
        return '%s object <%s>' % (self.__class__.__name__, self.name)
