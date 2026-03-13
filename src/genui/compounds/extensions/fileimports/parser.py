from abc import ABC, abstractmethod
from typing import Iterable


class FileParser(ABC):

    def __init__(self, file, molset):
        self.molset = molset
        self.file = file
        self.path = file.path

    @abstractmethod
    def parse(self) -> Iterable[tuple]:
        pass