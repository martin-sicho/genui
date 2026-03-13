import json
import logging
import traceback

import numpy as np
from django.core.files.base import ContentFile
from django.db import models, transaction
from django.db.models import Avg

from genui.compounds.models import MolSet, Molecule
from genui.models.models import Model, TrainingStrategy, ModelFile
from genui.qsar.models import EmbeddingCalculator

from django_rdkit import models as djrdkit

logger = logging.getLogger(__name__)


class Map(Model):
    molsets = models.ManyToManyField(MolSet, related_name="maps")

    @property
    def chemspaceJSON(self):
        query = self.files.filter(kind=ModelFile.AUXILIARY, note__exact='chemspaceJSON')
        if query:
            return query.get()
        else:
            return None


    def getChemSpaceJSDict(self, properties=("AMW", "NUMHEAVYATOMS", "NUMAROMATICRINGS", "HBA", "HBD", "LOGP", "TPSA",), file=None):
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
            f'{algorithm.name}-y',
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
                prop_calculator = getattr(djrdkit, prop, None)
                if prop_calculator:
                    molecule = molecule.annotate(**{lookup: prop_calculator('entity__rdMol')})
            molecule = molecule.first()
            if molecule:
                for prop in properties:
                    val = getattr(molecule, f"rdkit_prop_{prop}", None)
                    point_data['features'].append(val)

                # attach activities
                for activity_type in activity_types:
                    query = molecule.activities.filter(type__value=activity_type, source__in=activity_sets.values())
                    if query.exists():
                        try:
                            avg = query.aggregate(Avg('value'))
                            if avg is None or 'value__avg' not in avg:
                                raise ValueError(f"Invalid aggregate result for activity type {activity_type}")
                            avg_value = avg['value__avg']
                            point_data['features'].append(avg_value)
                        except Exception as e:
                            logger.error(f"Error calculating average for activity type {activity_type}: {str(e)}"
                                         f" {traceback.format_exc()}")
                            point_data['features'].append(None)
                    else:
                        point_data['features'].append(None)
                point_data['features'] = [p if p is not None and not np.isnan(p) else None for p in point_data['features']]
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
        try:
            if self.chemspaceJSON:
                self.chemspaceJSON.delete()
            model_file = ModelFile.create(
                self,
                f'chemspacejs.json',
                ContentFile(''),
                kind=ModelFile.AUXILIARY,
                note='chemspaceJSON',
            )
            chemspace_dict = self.getChemSpaceJSDict(file=model_file.path)
            if not chemspace_dict:
                raise ValueError("Failed to generate ChemSpace dictionary")
            print('Done.')
            return model_file
        except Exception as e:
            print(f"Error saving ChemSpace JSON for Map {self.pk}: {str(e)}")
            logger.error(f"Error saving ChemSpace JSON for Map {self.pk}: {str(e)}", exc_info=True)
            raise


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
    embeddings = models.ManyToManyField(EmbeddingCalculator)

