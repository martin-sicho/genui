from django.db import models
from polymorphic.models import PolymorphicModel
from abc import ABCMeta, abstractmethod

from django.utils import timezone

class PolymorphicAbstractModelMeta(ABCMeta, type(PolymorphicModel)):
    pass

class PolymorphicAbstractModel(PolymorphicModel):
    __metaclass__ = PolymorphicAbstractModelMeta

    class Meta:
        abstract = True

class BaseProject(PolymorphicAbstractModel):

    name = models.CharField(max_length=256, blank=False)
    description = models.TextField(max_length=10000, blank=True)
    created = models.DateTimeField(blank=True)
    updated = models.DateTimeField(blank=True, verbose_name="Last Update")

    class Meta:
        abstract = True

    @abstractmethod
    def update(self):
        pass

class BaseDataSet(PolymorphicAbstractModel):

    project = models.ForeignKey(BaseProject, on_delete=models.CASCADE)
    name = models.CharField(max_length=256, blank=False)
    description = models.TextField(max_length=10000, blank=True)
    created = models.DateTimeField()
    updated = models.DateTimeField('last_updated')

    class Meta:
        abstract = True

    @abstractmethod
    def update(self):
        pass

class Project(BaseProject):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.pk is None:
            self.created = timezone.now()
            self.update()

    def update(self):
        self.updated = timezone.now()

    def save(self, *args, **kwargs):
        self.update()
        super().save(*args, **kwargs)

class DataSet(BaseDataSet):

    project = models.ForeignKey(Project, on_delete=models.CASCADE, blank=False)

    class Meta:
        abstract = True

    def __str__(self):
        return '%s object (%s)' % (self.__class__.__name__, self.name)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.pk is None:
            self.created = timezone.now()
            self.update()

    def update(self):
        if self.pk is not None:
            self.project.update()
        self.updated = timezone.now()

    def save(self, *args, **kwargs):
        self.update()
        self.project.save()
        super().save(*args, **kwargs)
