from abc import ABC, abstractmethod

from .controls import BoundaryControl
from .initial_data import InitialData
from .params import StringParams

class ControlDesigner(ABC):
    name = "построитель"
    description = ""
    options: dict = {}

    @abstractmethod
    def build(self, params: StringParams, data: InitialData, horizon: float, **options) -> BoundaryControl:
        """Возвращает управление на (0, horizon)."""

REGISTRY: dict = {}

def register(designer: ControlDesigner):
    REGISTRY[designer.name] = designer
    return designer
