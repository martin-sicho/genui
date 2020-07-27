"""
chembl

Created by: Martin Sicho
On: 18-12-19, 14:38
"""
import traceback
import json

from django.db import transaction, IntegrityError

from genui.compounds.initializers.exceptions import SMILESParsingError
from genui.compounds.initializers.base import MolSetInitializer
from genui.compounds.models import ActivityUnits, ActivityTypes
from genui.utils.exceptions import GenUIException
from . import models


class ChEMBLSetInitializer(MolSetInitializer):

    def __init__(self, instance: models.ChEMBLCompounds, progress_recorder=None, targets=tuple(), max_per_target=None):
        super().__init__(instance, progress_recorder=progress_recorder)
        from chembl_webresource_client.new_client import new_client
        self.CHEMBL_ACTIVITIES = new_client.activity
        self.extracted_fields=(
            "MOLECULE_CHEMBL_ID"
            , "PARENT_MOLECULE_CHEMBL_ID"
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

    def createMolecule(self, entity, molecule_class, create_kwargs=None):
        chemblID = create_kwargs['chemblID']
        if molecule_class.objects.filter(chemblID=chemblID).exists():
            return molecule_class.objects.get(chemblID=chemblID)
        else:
            return super().createMolecule(entity, molecule_class, create_kwargs)

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
        current_item = 0
        for target, query in zip(self.targets, queries):
            target = models.ChEMBLTarget.objects.get_or_create(targetID=target)[0]
            print(f'Processing target: {target.targetID}')
            current_target_count += 1
            for result in query:
                current_item += 1
                if not self.max_per_target and self.progress_recorder:
                    self.progress_recorder.set_progress(current_item, progress_total)
                elif self.progress_recorder:
                    self.progress_recorder.set_progress(self.unique_mols, progress_total)

                # move on if we reached the maximum number of molecules per target in the set
                if self.max_per_target and (self.unique_mols / current_target_count) >= self.max_per_target:
                    break

                try:
                    compound_data = dict()
                    for field in self.extracted_fields:
                        compound_data[field] = result[field.lower()]
                    if not compound_data["CANONICAL_SMILES"]:
                        raise SMILESParsingError("", None, "Missing canonical SMILES string for molecule: {0}".format(compound_data["MOLECULE_CHEMBL_ID"]))
                except SMILESParsingError as exp:
                    traceback.print_exc()
                    self.errors.append(exp)
                    continue

                # check if this molecule is the parent molecule, if not skip it
                mol_chembl_id = compound_data['MOLECULE_CHEMBL_ID']
                parent_chembl_id = compound_data['MOLECULE_CHEMBL_ID']
                if mol_chembl_id != parent_chembl_id:
                    continue

                # create the molecule object and attach it to the instance
                try:
                    print(f"Creating {mol_chembl_id}...")
                    molecule = self.addMoleculeFromSMILES(compound_data["CANONICAL_SMILES"], models.ChEMBLMolecule, {"chemblID" : mol_chembl_id})
                    print(f"{mol_chembl_id} saved. Progress: {self.unique_mols}/{progress_total}")
                except GenUIException as exp:
                    print(f"The following exception happened while processing {mol_chembl_id}: {json.dumps(exp.asJSON(), indent=4)}")
                    traceback.print_exc()
                    self.errors.append(exp)
                    continue
                except IntegrityError as exp:
                    print(f"Database Integrity violation while creating molecule: {mol_chembl_id}")
                    traceback.print_exc()
                    self.errors.append(GenUIException(exp, f"Database Integrity violation while creating molecule: {mol_chembl_id}"))
                    continue

                # add found assay into assays or skip unwanted assays
                assay = models.ChEMBLAssay.objects.get_or_create(assayID=compound_data['ASSAY_CHEMBL_ID'])[0]
                with transaction.atomic():
                    # check if there are activity data
                    if compound_data['STANDARD_VALUE'] is None:
                        self.errors.append(GenUIException(None, f'No standard activity value found for molecule "{mol_chembl_id}" in assay "{assay.assayID}"')) # TODO: make a specific exception
                        continue

                    # add activity data
                    type_ = None
                    if compound_data['STANDARD_VALUE']:
                        units_val = compound_data['STANDARD_UNITS'] if compound_data['STANDARD_UNITS'] else 'Unknown'
                        type_val = compound_data['STANDARD_TYPE'] if compound_data['STANDARD_TYPE'] else 'Unknown'
                        units = ActivityUnits.objects.get_or_create(value=units_val)[0]
                        type_ = ActivityTypes.objects.get_or_create(value=type_val)[0]
                        activity = models.ChEMBLActivity.objects.create(
                            value=compound_data['STANDARD_VALUE'],
                            units=units,
                            source=self.activities,
                            molecule=molecule,
                            type=type_,
                            relation=compound_data['STANDARD_RELATION'] if compound_data['STANDARD_RELATION'] else "= (auto-assigned)",
                            assay=assay,
                            target=target,
                            comment=compound_data['ACTIVITY_COMMENT']
                        )
                        activity.save()
                    if type_ and activity and compound_data['PCHEMBL_VALUE']:
                        pchembl_value = models.ChEMBLActivity.objects.create(
                            value=compound_data['PCHEMBL_VALUE'],
                            source=self.activities,
                            molecule=molecule,
                            type=ActivityTypes.objects.get_or_create(
                                value=f'{type_.value}_pChEMBL'
                            )[0],
                            relation=compound_data['STANDARD_RELATION'] if compound_data['STANDARD_RELATION'] else "= (auto-assigned)",
                            assay=assay,
                            target=target,
                            comment=compound_data['ACTIVITY_COMMENT'],
                            parent=activity
                        )
                        pchembl_value.save()
        return self.unique_mols

    def updateInstance(self):
        self.instance.activities.all().delete()
        self.instance.molecules.clear()
        self.populateInstance()