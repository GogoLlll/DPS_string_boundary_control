from dataclasses import dataclass

@dataclass(frozen=True)
class StringParams:

    l: float
    a: float

    def __post_init__(self):
        if not (self.l > 0):
            raise ValueError("Длина струны l должна быть положительной")
        if not (self.a > 0):
            raise ValueError("Скорость волны a должна быть положительной")

    @property
    def t0(self) -> float:
        return self.l / self.a
