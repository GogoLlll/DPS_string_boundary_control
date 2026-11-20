from dataclasses import dataclass, field

from .controls import BoundaryControl, FunctionControl
from .functions import Expr
from .initial_data import InitialData
from .params import StringParams
from .solution import Solution

@dataclass
class Problem:
    params: StringParams
    data: InitialData
    control: BoundaryControl
    horizon: float
    texts: dict = field(default_factory=dict)
    _sources: object = field(default=None, repr=False, compare=False)

    def __post_init__(self):
        if not (self.horizon > 0):
            raise ValueError("Время T должно быть положительным")

    def solve(self) -> Solution:
        return Solution(self.params, self.data, self.control)

    def sources(self):
        """Источники скачков решения на (0, horizon); вычисляются один раз."""
        if self._sources is None:
            from .characteristics import find_sources
            self._sources = find_sources(self.params, self.data, self.control, self.horizon)
        return self._sources

    @classmethod
    def from_text(cls, l, a, T, y0, y1, u, active_until=None):
        """Собирает задачу из формул-строк. В формулах доступны параметры l, a, t0, T."""
        params = StringParams(float(l), float(a))
        T = float(T)
        constants = {"l": params.l, "a": params.a, "t0": params.t0, "T": T}
        y0_f, y1_f, u_f = Expr(y0, "x", constants), Expr(y1, "x", constants), Expr(u, "t", constants)
        data = InitialData.build(params, y0_f, y1_f)
        control = FunctionControl(u_f, active_until=active_until, label=str(u))
        return cls(params, data, control, T, texts={"y0": str(y0), "y1": str(y1), "u": str(u)})
