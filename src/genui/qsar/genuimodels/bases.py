from abc import ABC
from genui.utils.inspection import findSubclassByID, importFromPackage
from genui.qsar import models
from genui.utils.inspection import get_default_params


class EmbeddingCalculator(ABC):
    name = None
    no_init = True
    module = None
    _model = models.EmbeddingCalculator

    def __init__(self, builder):
        self.builder = builder
        self.instance = None

    def __call__(self, **kwargs):
        self.instance = getattr(self.module, self.name)(**kwargs)
        return self.instance

    def get_default_parameters(self):
        params =  get_default_params(self.name, self.module.__name__)
        if "kwargs" in params:
            params.pop("kwargs")
        return params

    @classmethod
    def getDjangoModel(cls, corePackage, update=False):
        if not cls.name:
            raise Exception('You have to specify a name for the embedding group in its class "name" property')

        ret, ret_created = cls._model.objects.get_or_create(name=cls.name)

        # just return if we are not setting up a new instance
        if not ret_created and not update:
            return ret

        if corePackage:
            ret.corePackage = corePackage
            ret.save()

        return ret


class EmbeddingBuilderMixIn:

    @staticmethod
    def findEmbeddingClass(name, corePackage="genui.qsar.genuimodels", subtype="embeddings"):
        module = importFromPackage(corePackage, subtype)
        return findSubclassByID(EmbeddingCalculator, module, "name", name)

    def __init__(self, instance: models.Model, progress=None, onFitCall=None):
        super().__init__(instance, progress, onFitCall)
        self.molsets = [self.instance.molset] if hasattr(self.instance, "molset") else self.instance.molsets.all()
        self.embeddingCalculators = [self._init_embedding_calculator(x) for x in self.training.embeddings.all()]

    def _init_embedding_calculator(self, django_model, subtype="embeddings"):
        class_ = self.findEmbeddingClass(django_model.name, django_model.corePackage, subtype)
        if "descriptors" in django_model.arguments:
            arguments = {class_.descriptors_name: django_model.arguments["descriptors"]}
            return class_(self)(**arguments) if class_ else None
        return class_(self)(**django_model.arguments) if class_ else None
