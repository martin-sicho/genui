import json

from django.core.files.base import ContentFile
from django.db import models, transaction
from django.db.models import Avg

from genui.compounds.models import MolSet, Molecule
from genui.models.models import Model, TrainingStrategy, ModelFile
from genui.qsar.models import DescriptorGroup

from django_rdkit import models as djrdkit


class Map(Model):
    molsets = models.ManyToManyField(MolSet, related_name="maps")

    @property
    def chemspaceJSON(self):
        query = self.files.filter(kind=ModelFile.AUXILIARY, note__exact='chemspaceJSON')
        if query:
            return query.get()
        else:
            return None


    def getChemSpaceJSDict(self, properties=("AMW",	"NUMHEAVYATOMS", "NUMAROMATICRINGS", "HBA", "HBD", "LOGP","TPSA",), file=None):
        ret = {
            'points': {},
            'compounds': {},
            'feature_names': [],
            'categories': [],
        }
        activity_sets = {}
        activity_types = set()
        molset_to_category = {}
        counter = 0
        for molset in self.molsets.all():
            category = {
                'id': molset.id,
                'label': molset.name,
                'points': []
            }
            ret['categories'].append(category)
            molset_to_category[molset.id] = counter
            counter += 1

            for activity_set in molset.activities.all():
                activity_sets[activity_set.id] = activity_set
                activity_types.update(activity_set.getActivityTypes())
        activity_types = sorted(list(activity_types))
        algorithm = self.trainingStrategy.algorithm
        ret['feature_names'] = [
            f'{algorithm.name}-x',
            f'{algorithm.name}-x',
        ] + list(properties) + activity_types


        for idx, point in enumerate(self.points.all()):
            point_data = {
                'features': [
                    point.x,
                    point.y,
                ]
            }
            molecule = Molecule.objects.filter(pk=point.molecule.id)

            # attach properties
            for prop in properties:
                lookup = f"rdkit_prop_{prop}"
                prop_calculator = getattr(djrdkit, prop)
                molecule = molecule.annotate(**{ lookup: prop_calculator('entity__rdMol')})
            molecule = molecule.get()
            for prop in properties:
                val = getattr(molecule, f"rdkit_prop_{prop}")
                point_data['features'].append(val)

            # attach activities
            for activity_type in activity_types:
                query = molecule.activities.filter(type__value=activity_type, source__in=activity_sets)
                if query:
                    avg = query.aggregate(Avg('value'))
                    point_data['features'].append(avg['value__avg'])
                else:
                    point_data['features'].append(None)

            ret['points'][point.id] = point_data
            ret['compounds'][point.id] = {
                'smiles': molecule.smiles,
                'id': molecule.id
            }
            for molset in molecule.providers.filter(pk__in=self.molsets.all()):
                ret['categories'][molset_to_category[molset.id]]['points'].append(point.id)

        if file:
            with open(file, mode='w', encoding='utf-8') as jsonfile:
                json.dump(ret, jsonfile)

        return ret

    @transaction.atomic
    def saveChemSpaceJSON(self):
        print(f'Saving ChemSpace.js JSON for {self.name}...')
        if self.chemspaceJSON:
            self.chemspaceJSON.delete()
        model_file = ModelFile.create(
            self,
            f'chemspacejs.json',
            ContentFile(''),
            kind=ModelFile.AUXILIARY,
            note='chemspaceJSON',
        )
        self.getChemSpaceJSDict(file=model_file.path)
        print('Done.')

        return model_file


class Point(models.Model):
    map = models.ForeignKey(Map, on_delete=models.CASCADE, null=False, related_name='points')
    molecule = models.ForeignKey(Molecule, on_delete=models.CASCADE, null=False)
    x = models.FloatField(blank=False, null=False)
    y = models.FloatField(blank=False, null=False)

    @property
    def smiles(self):
        return self.molecule.smiles

    @property
    def compoundSets(self):
        return self.molecule.providers.filter(id__in=[x.id for x in self.map.molsets.all()])

class MappingStrategy(TrainingStrategy):
    descriptors = models.ManyToManyField(DescriptorGroup)

