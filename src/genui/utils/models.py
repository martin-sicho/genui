"""
models

Created by: Martin Sicho
On: 1/12/20, 3:16 PM
"""
import os

from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.db import models


class OverwriteStorage(FileSystemStorage):

    def get_available_name(self, name, **kwargs):
        """Returns a filename that's free on the target storage system, and
        available for new content to be written to.

        Found at http://djangosnippets.org/snippets/976/

        This file storage solves overwrite on upload problem. Another
        proposed solution was to override the save method on the model
        like so (from https://code.djangoproject.com/ticket/11663):

        def save(self, *args, **kwargs):
            try:
                this = MyModelName.objects.get(id=self.id)
                if this.MyImageFieldName != self.MyImageFieldName:
                    this.MyImageFieldName.delete()
            except: pass
            super(MyModelName, self).save(*args, **kwargs)
        """
        # If the filename already exists, remove it as if it was a true file system
        if self.exists(name):
            os.remove(os.path.join(settings.MEDIA_ROOT, name))
        return name


def NON_POLYMORPHIC_CASCADE(collector, field, sub_objs, using):
    """
    This is a special cascade implementation to fix some delete errors
    when cascading polymorphic models.

    See: https://github.com/django-polymorphic/django-polymorphic/issues/229#issuecomment-398434412

    :param collector:
    :param field:
    :param sub_objs:
    :param using:
    :return:
    """

    return models.CASCADE(collector, field, sub_objs.non_polymorphic(), using)