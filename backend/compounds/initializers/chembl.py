"""
chembl

Created by: Martin Sicho
On: 18-12-19, 14:38
"""
import time
import traceback

from tqdm import tqdm
from django.db import transaction

from compounds.initializers.exceptions import SMILESParsingError
from .base import MolSetInitializer
from ..models import ChEMBLCompounds, ChEMBLMolecule, ChEMBLAssay, ChEMBLActivities, ChEMBLActivity, ActivityUnits, \
    ChEMBLTarget
from chembl_webresource_client.new_client import new_client

class ChEMBLSetInitializer(MolSetInitializer):

    def __init__(self, instance: ChEMBLCompounds, progress_recorder=None, targets=tuple(), max_per_target=None):
        super().__init__(instance, progress_recorder=progress_recorder)
        self.CHEMBL_ACTIVITIES = new_client.activity
        self.extracted_fields=(
            "MOLECULE_CHEMBL_ID"
            , "CANONICAL_SMILES"
            , 'ASSAY_CHEMBL_ID'
            , 'TARGET_CHEMBL_ID'
            , "PCHEMBL_VALUE"
            , "ACTIVITY_COMMENT"
            , "STANDARD_UNITS"
            , "STANDARD_TYPE"
            , "STANDARD_RELATION"
            , "STANDARD_VALUE"
        )
        self.activities = None
        activities = list(instance.activities.all())
        if not activities:
            self.activities = ChEMBLActivities.objects.create(name=f"{instance.name} - activities", description="Auto-assigned set of activities imported from ChEMBL.", project=instance.project, molecules=instance)
        else:
            self.activities = activities

        if targets:
            self.targets = targets
        else:
            self.targets = [x.targetID for x in instance.targets.all()]

        self.errors = []
        self.max_per_target = max_per_target

    def populateInstance(self):
        queries = []
        for target in self.targets:
            query = self.CHEMBL_ACTIVITIES.filter(
                target_chembl_id=target
                , target_organism='Homo sapiens'
                # , pchembl_value__isnull=False
            )
            queries.append(query)
        # queries = [x[0:self.max_per_target] if self.max_per_target else x for x in queries]
        progress_total = sum(len(x) for x in queries) if not self.max_per_target else self.max_per_target * len(self.targets)
        current_target_count = 0
        for target, query in zip(self.targets, queries):
            target = ChEMBLTarget.objects.get_or_create(targetID=target)[0]
            print(f'Processing target: {target.targetID}')
            current_target_count += 1
            for result in query:

                # move on if we reached the maximum number of molecules per target in the set
                if self.max_per_target and (self.unique_mols / current_target_count) >= self.max_per_target:
                    break

                try:
                    compound_data = dict()
                    for field in self.extracted_fields:
                        compound_data[field] = result[field.lower()]
                    if not compound_data["CANONICAL_SMILES"]:
                        raise SMILESParsingError("", "Missing canonical SMILES string for molecule: {0}".format(compound_data["MOLECULE_CHEMBL_ID"]))
                except SMILESParsingError as exp:
                    traceback.print_exc()
                    self.errors.append(exp)
                    continue

                # create the molecule object
                mol_chembl_id = compound_data['MOLECULE_CHEMBL_ID']
                with transaction.atomic():
                    molecule = self.addMoleculeFromSMILES(compound_data["CANONICAL_SMILES"], ChEMBLMolecule, {"chemblID" : mol_chembl_id})
                    print(f"Saving {mol_chembl_id}... ({self.unique_mols}/{progress_total})")
                    molecule.save()

                # add found assay into assays or skip unwanted assays
                assay = ChEMBLAssay.objects.get_or_create(assayID=compound_data['ASSAY_CHEMBL_ID'])[0]
                with transaction.atomic():
                    # add activity data
                    if compound_data['STANDARD_VALUE'] and compound_data['STANDARD_RELATION'] and compound_data['STANDARD_TYPE']:
                        activity = ChEMBLActivity.objects.create(
                            value=compound_data['STANDARD_VALUE'],
                            units=ActivityUnits.objects.get_or_create(value=compound_data['STANDARD_UNITS'])[0] if compound_data['STANDARD_UNITS'] else None,
                            source=self.activities,
                            molecule=molecule,
                            type=compound_data['STANDARD_TYPE'],
                            relation=compound_data['STANDARD_RELATION'],
                            assay=assay,
                            target=target,
                            comment=compound_data['ACTIVITY_COMMENT']
                        )
                        activity.save()
                    if compound_data['PCHEMBL_VALUE']:
                        pchembl_value = ChEMBLActivity.objects.create(
                            value=compound_data['PCHEMBL_VALUE'],
                            source=self.activities,
                            molecule=molecule,
                            type="PCHEMBL_VALUE",
                            relation=compound_data['STANDARD_RELATION'],
                            assay=assay,
                            target=target,
                            comment=compound_data['ACTIVITY_COMMENT']
                        )
                        pchembl_value.save()
                if self.progress_recorder:
                    self.progress_recorder.set_progress(self.unique_mols, progress_total)
        return self.unique_mols

    def updateInstance(self):
        if True:
            total = 60
            for i in range(total):
                print(i)
                time.sleep(1)
                if self.progress_recorder:
                    self.progress_recorder.set_progress(i, total)
            return total
        else:
            # FIXME: make this happen
            instance = self.getInstance()
            instance.activities.clear()
            instance.molecules.clear()
            self.populateInstance()