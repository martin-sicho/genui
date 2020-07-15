"""
initializer

Created by: Martin Sicho
On: 7/13/20, 4:12 PM
"""

from genui.compounds.models import Activity, ActivityTypes, ActivityUnits, ActivitySet
from genui.compounds.extensions.fileimports.initializer import FileInitializer
from .models import SDFMolecule


class SDFSetInitializer(FileInitializer):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.activitySet = None

    def parserCallback(self, smile, props):
        molecule = self.addMoleculeFromSMILES(smile, SDFMolecule, {'name' : props['name']})
        print(f'Imported molecule: {molecule.name}')

        # attach activities
        if self.instance.activitiesProp in props and self.instance.activityTypesProp in props:
            print(f"Found activities for molecule: {molecule.name}")
            if not self.activitySet:
                self.activitySet = ActivitySet.objects.create(
                    name=f'{self.instance.name} (activities from file)',
                    project=self.instance.project,
                    molecules=self.instance,
                )

            sep = self.instance.dataSeparator
            values = [float(x) for x in props[self.instance.activitiesProp].split(sep)]
            types = props[self.instance.activityTypesProp].split(sep)
            assert len(values) == len(types) # TODO: this should be a custom exception
            units = props[self.instance.activityUnitsProp].split(sep) if self.instance.activityUnitsProp in props else len(types) * [None]

            for value, _type, unit in zip(values, types, units):
                unit = None if not unit or unit == 'None' else ActivityUnits.objects.get_or_create(value=unit)[0]
                _type = ActivityTypes.objects.get_or_create(value=_type)[0]
                Activity.objects.create(
                    value=value,
                    type=_type,
                    units=unit,
                    source=self.activitySet,
                    molecule=molecule
                )
                print(f"Imported activity value: {_type.value} = {value} {unit.value if unit else ''}")


    def populateInstance(self):
        self.parseMols(self.parserCallback)
        return self.unique_mols

