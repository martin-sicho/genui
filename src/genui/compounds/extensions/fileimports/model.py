from genui.compounds.models import MolSet


class FileCompounds(MolSet):

    class Meta:
        abstract = True

    @property
    def file(self):
        return self.files.all()[0].file