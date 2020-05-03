from django.db import models

# Create your models here.
from djcelery_model.models import TaskMixin

from genui.utils.extensions.tasks.models import TaskShortcutsMixIn, PolymorphicTaskManager
from genui.compounds.models import MolSet
from genui.projects.models import DataSet


class Generator(TaskShortcutsMixIn, TaskMixin, DataSet):
    objects = PolymorphicTaskManager()

    # TODO: it would be useful to have this as an abstract method if possible
    def get(self, n_samples) -> [str]:
        """
        All generators should override this. This should
        return a list of SMILES strings of size N.

        :return: list of SMILES strings
        """

        raise NotImplemented("You have to override this method in subclasses.")

class GeneratedMolSet(MolSet):
    source = models.ForeignKey(Generator, on_delete=models.CASCADE, null=False, related_name="compounds")
