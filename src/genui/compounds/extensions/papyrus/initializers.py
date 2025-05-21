import traceback
import json
import math

from django.db import transaction, IntegrityError

from genui.compounds.initializers.exceptions import SMILESParsingError, CompoundImportException
from genui.compounds.initializers.base import MolSetInitializer
from genui.compounds.models import ActivityUnits, ActivityTypes
from . import models

from papyrus_scripts.reader import read_papyrus, read_protein_set
from papyrus_scripts.preprocess import (consume_chunks, keep_organism)

class PapyrusSetInitializer(MolSetInitializer):

    def __init__(self, instance: models.PapyrusCompounds, progress_recorder=None, targets=tuple(), max_per_target=None):
        super().__init__(instance, progress_recorder=progress_recorder)
        self.CHUNK_SIZE = 30_000
        self.total_chunks = math.ceil(707_461/self.CHUNK_SIZE) # 707461 compound-protein activity points without stereochemistry in Papyrus++, 59775912 in full dataset
        self.compound_data = read_papyrus(is3d=False, chunksize=self.CHUNK_SIZE, source_path=None, plusplus=True, version='latest')
        self.target_data = read_protein_set(source_path=None) #for additional target information, add later
        self.activities = self.initActivities()
        self.errors = []
        self.max_per_target = max_per_target

        # FIXME no targets entered on first run (fix on GUI level)
        if targets:
            self.targetIDs = targets
        else:
            self.targetIDs = [x.targetID for x in instance.targets.all()]

    def initActivities(self):
        activities = list(self.instance.activities.all())
        if not activities:
            activities = models.PapyrusActivities.objects.create(name=f"{self.instance.name} Activities (imported)", description=f"Activity information downloaded from ChEMBL when the {self.instance.name} compound set was created.", project=self.instance.project, molecules=self.instance)
        return activities
 
    def createMolecule(self, entity, molecule_class, create_kwargs=None):
        # FIXME CID not unique in dataset, some molecules could be skipped and liked incorrectly
        CID = create_kwargs['CID']
        if molecule_class.objects.filter(CID=CID).exists():
            return molecule_class.objects.get(CID=CID)
        else:
            return super().createMolecule(entity, molecule_class, create_kwargs)
        
    def populateInstance(self):
        target_counts = dict.fromkeys(self.targetIDs, 0)
        for i_chunk, chunk in enumerate(self.compound_data):
            chunk = keep_organism(data=chunk, protein_data=self.target_data, organism=['Human'], generic_regex=True)
            for targetID in self.targetIDs:
                target = models.PapyrusTarget.objects.get_or_create(targetID=targetID)[0]
                # move on if we reached the maximum number of molecules per target in the set
                if self.max_per_target and target_counts[targetID] >= self.max_per_target:
                    continue

                query = chunk[chunk['target_id'] == targetID+'_WT']

                for i_query, row in enumerate(query.itertuples(index=False)):
                    # Excluding non wild type targets
                    if row.Protein_Type != 'WT':
                        print('Log non WT, id: ', row.Activity_ID)
                        continue
                    # move on if we reached the maximum number of molecules per target in the set
                    if self.max_per_target and target_counts[targetID] >= self.max_per_target:
                        break
                    target_counts[targetID] = target_counts[targetID] + 1

                    if self.progress_recorder:
                        if self.max_per_target:
                            progress_total = self.max_per_target * len(self.targetIDs)
                            self.progress_recorder.set_progress(sum(target_counts.values()), progress_total)
                        else:
                            chunk_progress = (i_chunk/self.total_chunks)*94
                            molecule_progress = (i_query/query.shape[0])*6
                            progress = int(chunk_progress+molecule_progress)
                            self.progress_recorder.set_progress(progress, 100)

                    try:
                        if not row.SMILES:
                            raise SMILESParsingError("", None, "Missing SMILES string for molecule: {0}".format(row.CID.split(";")[0]))
                    except SMILESParsingError as exp:
                        traceback.print_exc()
                        self.errors.append(exp)
                        continue
                    
                    # create the molecule object and attach it to the instance
                    CID = row.CID
                    try:
                        print(f"Creating {CID}...")
                        molecule = self.addMoleculeFromSMILES(row.SMILES, models.PapyrusMolecule, {"CID" : CID})
                        print(f"{CID} saved.") # more useful message
                    except CompoundImportException as exp:
                        print(f"The following exception happened while processing {CID}: {json.dumps(exp.asJSON(), indent=4)}")
                        traceback.print_exc()
                        self.errors.append(exp)
                        continue
                    except IntegrityError as exp:
                        print(f"Database Integrity violation while creating molecule: {CID}")
                        traceback.print_exc()
                        self.errors.append(CompoundImportException(exp, f"Database Integrity violation while creating molecule: {CID}"))
                        continue

                    # add found assay into assays or skip unwanted assays
                    # FIXME multiple assays
                    assay = models.PapyrusAssay.objects.get_or_create(assayID=row.AID)[0]
                    with transaction.atomic():
                        # check if there are activity data
                        if row.pchembl_value is None:
                            self.errors.append(CompoundImportException(None, f'No activity value found for molecule "{CID}" in assay "{assay.assayID}"')) # TODO: make a specific exception
                            continue

                        # add activity data
                        type_ = None
                        if row.pchembl_value:
                            ic = list(map(int, row.type_IC50.split(";")))
                            ec = list(map(int, row.type_EC50.split(";")))
                            kd = list(map(int, row.type_KD.split(";")))
                            ki = list(map(int, row.type_Ki.split(";")))
                            values = list(map(float, row.pchembl_value.split(";")))
                            transposed = [list(x) for x in zip(ic, ec, kd, ki)]
                            type_value_pairs = []
                            for i, x in enumerate(transposed):
                                type_val =  'IC50' if x[0] == 1 else \
                                            'EC50' if x[1] == 1 else \
                                            'KD' if x[2] == 1 else \
                                            'Ki' if x[3] == 1 else \
                                            'Unknown'
                                type_value_pairs.append((type_val, values[i]))
                            #units = ActivityUnits.objects.get_or_create(value='Unknown')[0]
                            for pair in type_value_pairs:
                                type_ = ActivityTypes.objects.get_or_create(value=pair[0])[0]
                                value = pair[1]
                                activity = models.PapyrusActivity.objects.create(
                                    value = value,
                                    #units = units,
                                    source = self.activities,
                                    molecule = molecule,
                                    type = type_,
                                    relation = row.relation if row.relation else "= (auto-assigned)",
                                    assay = assay,
                                    target = target,
                                )
                                activity.save()
                            # add mean and median of pchembl_value (recorded in Papyrus)
                            mean = models.PapyrusActivity.objects.create(
                                    value = row.pchembl_value_Mean,
                                    source = self.activities,
                                    molecule = molecule,
                                    type = ActivityTypes.objects.get_or_create(value='Mean')[0],
                                    relation = '=',
                                    assay = assay,
                                    target = target,
                                )
                            mean.save()
                            median = models.PapyrusActivity.objects.create(
                                    value = row.pchembl_value_Median,
                                    source = self.activities,
                                    molecule = molecule,
                                    type = ActivityTypes.objects.get_or_create(value='Median')[0],
                                    relation = '=',
                                    assay = assay,
                                    target = target,
                                )
                            median.save()
        return self.unique_mols

    def updateInstance(self):
        self.progress_recorder.set_progress(0, 100, description='Deleting existing records.')
        self.instance.activities.all().delete()
        self.activities = self.initActivities()
        self.instance.molecules.clear()
        self.populateInstance()