import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from wavesys.problem import Problem  # noqa: E402

@pytest.fixture
def smooth():
    return Problem.from_text(
        l=1.5, a=2.0, T=3.0,
        y0="0.3*(1 - x/l) + 0.7*sin(pi*x/l) + 0.2*sin(3*pi*x/l)",
        y1="0.5*sin(2*pi*x/l)",
        u="0.3*cos(3*t) + 0.4*sin(2*t)^2")

@pytest.fixture
def inner_points():
    return (np.arange(80) + 0.37) / 80
