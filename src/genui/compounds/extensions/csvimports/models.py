from django.db import models

from genui.compounds.extensions.fileimports.model import FileCompounds
from genui.compounds.models import Molecule


class CSVCompounds(FileCompounds):
    nameCol = models.CharField(max_length=256, null=True, default='NAME', blank=False)
    smilesCol = models.CharField(max_length=256, null=False, default='SMILES', blank=False)
    activityCol = models.CharField(max_length=256, null=False, default='ACTIVITY', blank=False)
    activityTypeCol = models.CharField(max_length=256, null=False, default='ACTIVITY_TYPE', blank=False)
    activityUnitsCol = models.CharField(max_length=256, null=True, default='ACTIVITY_UNITS', blank=False)
    colSeparator = models.CharField(max_length=1, null=False, default=',', blank=False)
    emptyValue = models.CharField(max_length=256, null=False, default='NA', blank=True)

class CSVMolecule(Molecule):
    name = models.CharField(max_length=1024, null=False, blank=True)

    def __str__(self):
        return '%s object <%s>' % (self.__class__.__name__, self.name)
