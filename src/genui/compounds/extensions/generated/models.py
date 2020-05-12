from django.db import models

# Create your models here.
from genui.compounds.models import MolSet
from genui.generators.models import Generator


class GeneratedMolSet(MolSet):
    source = models.ForeignKey(Generator, on_delete=models.CASCADE, null=False, related_name="compounds")