"""
initializer

Created by: Martin Sicho
On: 7/16/20, 1:54 PM
"""
from genui.compounds.extensions.fileimports.initializer import FileInitializer
from genui.compounds.models import ActivityUnits, ActivityTypes, Activity, ActivitySet
from genui.utils.exceptions import GenUIException
from .models import CSVMolecule

class CSVParsingError(GenUIException):
    pass

class CSVSetInitializer(FileInitializer):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.activitySet = None

    def parserCallback(self, smile, props):
        molecule = self.addMoleculeFromSMILES(smile, CSVMolecule, {'name' : props['name']})
        print(f'Imported molecule: {molecule.name}')

        # attach activities
        if 'activityType' in props and 'activity' in props:
            value = props['activity']
            if not value or value == self.instance.emptyValue or type(value) not in (int, float):
                self.errors.append(CSVParsingError(None, f'Failed to parse an activity value ({value}) for molecule {molecule.name} : {smile}'))
                return
            value = float(value)

            _type = props['activityType']
            if not _type or _type == self.instance.emptyValue or type(_type) != str:
                self.errors.append(CSVParsingError(None, f'Failed to parse activity type ({_type}) for molecule {molecule.name} : {smile}'))
                return
            _type = ActivityTypes.objects.get_or_create(value=_type)[0]

            if not self.activitySet:
                self.activitySet = ActivitySet.objects.create(
                    name=f'{self.instance.name} (activities from file)',
                    project=self.instance.project,
                    molecules=self.instance,
                )

            units = props['activityUnits']
            units = None if not units or type(units) != str else ActivityUnits.objects.get_or_create(value=units)[0]

            Activity.objects.create(
                value=value,
                type=_type,
                units=units,
                source=self.activitySet,
                molecule=molecule
            )
            print(f"Imported activity value: {_type.value} = {value} {units.value if units else ''}")

    def populateInstance(self):
        self.parseMols(self.parserCallback)
        return self.unique_mols

