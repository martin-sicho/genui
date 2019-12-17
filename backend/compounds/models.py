from django.db import models
from projects.models import BaseDataSet, Project

from django.utils import timezone

# Create your models here.

class MolSet(BaseDataSet):

    project = models.ForeignKey(Project, on_delete=models.CASCADE)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.pk is None:
            self.created = timezone.now()
            self.update()

    def update(self):
        self.project.update()
        self.updated = timezone.now()

    def save(self, *args, **kwargs):
        self.update()
        self.project.save()
        super().save(*args, **kwargs)
