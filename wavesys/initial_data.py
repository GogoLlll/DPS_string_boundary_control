from dataclasses import dataclass
from typing import Callable

import numpy as np

from .functions import Primitive
from .params import StringParams

@dataclass
class InitialData:
    y0: Callable
    y1: Callable
    F: Callable

    @classmethod
    def build(cls, params: StringParams, y0, y1, F=None):
        if F is None:
            F = Primitive(y1, 0.0, params.l)
        return cls(y0=y0, y1=y1, F=F)

    @classmethod
    def zero(cls):
        z = lambda x: np.zeros(np.shape(x))
        return cls(y0=z, y1=z, F=z)
