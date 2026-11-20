import numpy as np
import pytest

from wavesys import results
from wavesys.functions import ExpressionError
from wavesys.presets import DEFAULT, PRESETS
from wavesys.problem import Problem

def test_from_text_exposes_parameters_in_formulas():
    pr = Problem.from_text(l=2.0, a=4.0, T=3.0, y0="x/l", y1="a + t0", u="T*t")
    assert pr.params.t0 == 0.5
    assert pr.data.y0(1.0) == pytest.approx(0.5)
    assert pr.data.y1(0.3) == pytest.approx(4.5)
    assert pr.control(np.array([0.2]), 0.5)[0] == pytest.approx(0.6)

def test_numeric_primitive_matches_exact():
    pr = Problem.from_text(l=2.0, a=1.0, T=1.0, y0="0", y1="cos(2*x)", u="0")
    x = np.linspace(0, 2, 41)
    assert np.allclose(pr.data.F(x), np.sin(2 * x) / 2, atol=1e-13)

@pytest.mark.parametrize("kw", [dict(l=-1), dict(a=0), dict(T=0), dict(T=-2)])
def test_rejects_bad_numbers(kw):
    args = dict(l=1.0, a=1.0, T=1.0, y0="0", y1="0", u="0")
    args.update(kw)
    with pytest.raises(ValueError):
        Problem.from_text(**args)

def test_rejects_wrong_variable():
    with pytest.raises(ExpressionError):
        Problem.from_text(l=1, a=1, T=1, y0="sin(t)", y1="0", u="0")      # в y0 переменная - x
    with pytest.raises(ExpressionError):
        Problem.from_text(l=1, a=1, T=1, y0="0", y1="0", u="x")           # в u переменная - t

def test_sources_are_cached(smooth):
    assert smooth.sources() is smooth.sources()

def test_default_preset_exists():
    assert DEFAULT in PRESETS

@pytest.mark.parametrize("name", list(PRESETS))
def test_every_preset_solves_and_passes_checks(name):
    p = PRESETS[name]
    pr = Problem.from_text(**{k: p[k] for k in ("l", "a", "T", "y0", "y1", "u")})
    r = results.compute(pr, nx=41, nt=41, fd_nodes=60, norm_points=11)
    assert [c.name for c in r.checks if c.ok is False] == []
    assert np.all(np.isfinite(r.Y))

def test_left_moving_wave_preset_leaves_string_at_rest():
    p = PRESETS["Бегущая влево волна"]
    pr = Problem.from_text(**{k: p[k] for k in ("l", "a", "T", "y0", "y1", "u")})
    sol = pr.solve()
    assert sol.state_norm(0.0) > 0.01
    assert sol.state_norm(0.6 * pr.params.t0 + 1e-9) < 1e-12
