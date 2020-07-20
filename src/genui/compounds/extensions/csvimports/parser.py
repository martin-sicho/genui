"""
parser

Created by: Martin Sicho
On: 7/16/20, 2:17 PM
"""
from genui.compounds.extensions.fileimports.parser import FileParser
import pandas

class CSVParser(FileParser):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.counter = 0

    def processRow(self, row):
        self.counter += 1
        smiles = row[self.molset.smilesCol]
        activityUnits = row[self.molset.activityUnitsCol] if self.molset.activityUnitsCol in row.index else None
        name = row[self.molset.nameCol] if self.molset.nameCol in row.index else f'Molecule #{self.counter}'

        return smiles, {
            'activity' : row[self.molset.activityCol],
            'activityType' : row[self.molset.activityTypeCol],
            'activityUnits' : activityUnits,
            'name' : name,

        }

    def parse(self):
        df = pandas.read_csv(
            self.path,
            sep=self.molset.colSeparator,
            header=0,
            na_values=('', self.molset.emptyValue),
            keep_default_na=True
        )

        # TODO: this should be changed for something more efficient: https://stackoverflow.com/a/55557758
        return df.apply(self.processRow, axis=1).tolist()

