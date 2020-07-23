"""
exceptions

Created by: Martin Sicho
On: 7/23/20, 10:08 AM
"""
import json
import traceback


class GenUIException(Exception):

    def __init__(self, original, *args, **kwargs):
        super().__init__(*args)
        self.original = original

    def getData(self):
        return ''

    def __repr__(self):
        return self.asJSON()

    def asJSON(self):
        return json.dumps({
            "original" : str(type(self.original)) if self.original else '',
            "current" : str(type(self)),
            "reprOrig" : repr(self.original) if self.original else '',
            "tracebacks" : {
                "original" : traceback.extract_tb(self.original.__traceback__).format() if self.original else '',
                "current" : traceback.extract_tb(self.__traceback__).format()
            },
            "messages" : {
              "original" : [x for x in self.original.args] if self.original else [],
              "current" : [x for x in self.args]
            },
            "data" : self.getData()
        })


