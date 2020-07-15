"""
initializer

Created by: Martin Sicho
On: 7/13/20, 1:59 PM
"""
import traceback
from abc import abstractmethod

from genui.compounds.initializers.base import MolSetInitializer
from genui.utils.inspection import getObjectAndModuleFromFullName


class InvalidParserException(Exception):
    pass

class FileInitializer(MolSetInitializer):

    def __init__(self, instance, progress_recorder=None, parser_class = None, instance_file_attr = 'file'):
        super().__init__(instance, progress_recorder=progress_recorder)
        if not parser_class:
            raise InvalidParserException('No file parser specified for FileInitializer!')
        self.parserClass = getObjectAndModuleFromFullName(parser_class)[0]
        self.parser = self.parserClass(getattr(instance, instance_file_attr))

    def parseMols(self, callback):
        for item in self.parser.parse():
            smile = item[0]
            props = item[1]
            try:
                callback(smile, props)
            except Exception as exp:
                print(f"Callback error while adding molecule: {smile}")
                traceback.print_exc()
                self.errors.append(exp)

    @abstractmethod
    def populateInstance(self):
        pass

    def updateInstance(self):
        return self.populateInstance()