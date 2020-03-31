"""
chembl

Created by: Martin Sicho
On: 18-12-19, 14:38
"""
import time
import traceback

from django.db import transaction, IntegrityError

from compounds.initializers.exceptions import SMILESParsingError
from .base import MolSetInitializer
from .. import models

class ChEMBLSetInitializer(MolSetInitializer):

    def __init__(self, instance: models.ChEMBLCompounds, progress_recorder=None, targets=tuple(), max_per_target=None):
        super().__init__(instance, progress_recorder=progress_recorder)
        from chembl_webresource_client.new_client import new_client
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
            self.activities = models.ChEMBLActivities.objects.create(name=f"{instance.name} Activities (imported)", description=f"Activity information downloaded from ChEMBL when the {instance.name} compound set was created.", project=instance.project, molecules=instance)
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
            target = models.ChEMBLTarget.objects.get_or_create(targetID=target)[0]
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
                try:
                    molecule = self.addMoleculeFromSMILES(compound_data["CANONICAL_SMILES"], models.ChEMBLMolecule, {"chemblID" : mol_chembl_id})
                except IntegrityError as exp:
                    print(f"Integrity error for {mol_chembl_id}")
                    traceback.print_exc()
                    self.errors.append(exp)
                    continue
                print(f"Saving {mol_chembl_id}... ({self.unique_mols}/{progress_total})")
                molecule.save()

                # add found assay into assays or skip unwanted assays
                assay = models.ChEMBLAssay.objects.get_or_create(assayID=compound_data['ASSAY_CHEMBL_ID'])[0]
                with transaction.atomic():
                    # add activity data
                    units_val = compound_data['STANDARD_UNITS'] if compound_data['STANDARD_UNITS'] else 'Unknown'
                    type_val = compound_data['STANDARD_TYPE'] if compound_data['STANDARD_TYPE'] else 'Unknown'
                    if compound_data['STANDARD_VALUE']:
                        units = models.ActivityUnits.objects.get_or_create(value=units_val)[0]
                        type_ = models.ActivityTypes.objects.get_or_create(value=type_val)[0]
                        activity = models.ChEMBLActivity.objects.create(
                            value=compound_data['STANDARD_VALUE'],
                            units=units,
                            source=self.activities,
                            molecule=molecule,
                            type=type_,
                            relation=compound_data['STANDARD_RELATION'],
                            assay=assay,
                            target=target,
                            comment=compound_data['ACTIVITY_COMMENT']
                        )
                        activity.save()
                    if compound_data['PCHEMBL_VALUE']:
                        pchembl_value = models.ChEMBLActivity.objects.create(
                            value=compound_data['PCHEMBL_VALUE'],
                            source=self.activities,
                            molecule=molecule,
                            type=models.ActivityTypes.objects.get_or_create(
                                value='PCHEMBL'
                            )[0],
                            relation=compound_data['STANDARD_RELATION'],
                            assay=assay,
                            target=target,
                            comment=compound_data['ACTIVITY_COMMENT']
                        )
                        pchembl_value.save()
                    if not ((compound_data['STANDARD_VALUE'] is not None) | (compound_data['PCHEMBL_VALUE'] is not None)):
                        self.errors.append(Exception(f'No activity value found for molecule "{mol_chembl_id}" in assay "{assay.assayID}"')) # TODO: make a specific exception
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