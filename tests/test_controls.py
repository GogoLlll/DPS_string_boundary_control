import numpy as np

from wavesys.controls import FunctionControl, ZeroControl, find_jumps
from wavesys.functions import Expr


def test_zero_control():
    u = ZeroControl()
    assert np.all(u(np.linspace(0, 5, 11), 0.5) == 0)

def test_function_control_matches_function():
    u = FunctionControl(lambda t: 0.4 * t)
    t = np.linspace(0, 3, 301)
    assert np.allclose(u(t, 0.75), 0.4 * t, atol=1e-12)

def test_slot_is_shifted_time():
    u = FunctionControl(np.cos)
    s = np.linspace(0, 0.7, 8)
    assert np.allclose(u.slot(3, s, 0.75), np.cos(3 * 0.75 + s))

def test_active_until_switches_off_exactly():
    t0 = 0.08
    u = FunctionControl(lambda t: np.ones(np.shape(t)), active_until=4 * t0)
    assert np.all(u.slot(3, np.array([0.0, t0 * 0.999]), t0) == 1.0)
    assert np.all(u.slot(4, np.array([0.0, 0.01]), t0) == 0.0)        # в момент выключения уже ноль
    assert u.jump_times(t0, 1.0) == [4 * t0]

def test_find_jumps_ignores_kinks_and_steep_slopes():
    assert find_jumps(lambda t: np.abs(t - 0.4), 0, 1) == []
    assert find_jumps(lambda t: np.tanh(200 * (t - 0.5)), 0, 1) == []
    assert find_jumps(lambda t: np.sin(40 * t), 0, 1) == []

def test_find_jumps_locates_jumps():
    f = Expr("0.5*ind(t, 0.2, 0.55) + 0.1*t", "t")
    jumps = find_jumps(f, 0.0, 1.0)
    assert np.allclose(jumps, [0.2, 0.55], atol=1e-9)

def test_find_jumps_small_jump_on_large_background():
    jumps = find_jumps(lambda t: 100 * np.sin(t) + 1e-3 * (t > 0.7), 0.0, 1.0)
    assert np.allclose(jumps, [0.7], atol=1e-9)
