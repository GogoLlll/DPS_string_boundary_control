import numpy as np
import pytest

from wavesys.characteristics import (Source, find_sources, first_layer_region, fold, front_polyline,
                                     front_position, front_segments, point_info)
from wavesys.params import StringParams
from wavesys.problem import Problem

P = StringParams(2.0, 4.0)        # t0 = 0.5

@pytest.mark.parametrize("t_frac, x_frac, region", [
    (0.25, 0.5, "I"), (0.5, 0.25, "II"), (0.5, 0.75, "III"), (0.75, 0.5, "IV"),      # задача 2.3 лекций
])
def test_first_layer_regions(t_frac, x_frac, region):
    assert first_layer_region(P, t_frac * P.t0, x_frac * P.l) == region

def test_no_region_after_first_layer():
    assert first_layer_region(P, 1.2 * P.t0, 0.5 * P.l) is None

def test_point_info_segments_and_control_times():
    info = point_info(P, 2.6 * P.t0, 0.5 * P.l)          # z_A = z_B = 3.1 t0
    assert info["kA"] == 3 and info["kB"] == 3
    assert info["sA"] == pytest.approx(0.1 * P.t0)
    assert info["A_from"] == "-B0" and info["B_from"] == "-A0"
    assert np.allclose(info["control_times_A"], [1.1 * P.t0])                 # A3 = −B0 − 2 u1
    assert np.allclose(info["control_times_B"], [0.1 * P.t0, 2.1 * P.t0])     # B3 = −A0 + 2 u0 + 2 u2

def test_point_info_never_uses_future_control():
    rng = np.random.default_rng(1)
    for _ in range(200):
        t, x = rng.uniform(0, 6 * P.t0), rng.uniform(0, P.l)
        info = point_info(P, t, x)
        assert all(tau <= t + 1e-12 for tau in info["control_times_A"] + info["control_times_B"])

def test_fold_reflects_between_ends():
    assert np.allclose(fold(np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0]), 2.0), [0, 1, 2, 1, 0, 1])

def test_front_from_left_corner():
    src = Source("left", 0.0, 1.0)
    t = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    assert np.allclose(front_position(P, src, t), [0.0, 1.0, 2.0, 1.0, 0.0])

def test_front_from_right_corner_and_delayed_source():
    assert np.allclose(front_position(P, Source("right", 0.0, 1.0), np.array([0.0, 0.25, 0.5])), [2.0, 1.0, 0.0])
    pos = front_position(P, Source("left", 0.3, 1.0), np.array([0.1, 0.3, 0.55]))
    assert np.isnan(pos[0]) and np.allclose(pos[1:], [0.0, 1.0])

def test_front_polyline_has_vertices_at_reflections():
    t, x = front_polyline(P, Source("left", 0.0, 1.0), 1.3)
    assert np.allclose(t, [0.0, 0.5, 1.0, 1.3])
    assert np.allclose(x, [0.0, 2.0, 0.0, 1.2])

def test_sources_for_consistent_smooth_problem(smooth):
    assert find_sources(smooth.params, smooth.data, smooth.control, smooth.horizon) == []

def test_sources_detect_corners_and_control_jumps():
    pr = Problem.from_text(l=1.0, a=1.0, T=2.0, y0="0.5*(1 - x/l) + 0.2*x/l", y1="0", u="0.3*ind(t, 0, 0.4)")
    found = {(s.kind, round(s.time, 6)): s.size for s in pr.sources()}
    assert found[("left", 0.0)] == pytest.approx(0.3 - 0.5)      # u(0) − y0(0)
    assert found[("right", 0.0)] == pytest.approx(-0.2)          # −y0(l)
    assert found[("left", 0.4)] == pytest.approx(-0.3)           # скачок управления

def test_jump_of_solution_across_front_equals_source_size():
    """Скачок решения на прямом фронте равен δ, на отражённом — −δ (закон (XL) из reshenie.pdf)."""
    pr = Problem.from_text(l=1.0, a=1.0, T=2.0, y0="0.5*(1 - x/l)", y1="0", u="0.2")
    sol, (src,) = pr.solve(), pr.sources()
    x, e = 0.37, 1e-9
    assert sol.y(x + e, x) - sol.y(x - e, x) == pytest.approx(src.size, abs=1e-7)
    assert sol.y(2 - x + e, x) - sol.y(2 - x - e, x) == pytest.approx(-src.size, abs=1e-7)

def test_front_segments_merge_and_cancel():
    """Скачки через 2 t0 идут по одной линии и складываются; при нулевой сумме линия обрывается."""
    lines = front_segments(P, [Source("left", 0.0, 0.4), Source("left", 1.0, -0.4)], 3.0)      # 2 t0 = 1.0
    assert len(lines) == 1
    t, x, size = lines[0]
    assert np.allclose(t, [0.0, 0.5, 1.0]) and np.allclose(x, [0.0, 2.0, 0.0]) and size == pytest.approx(0.4)
    lines = front_segments(P, [Source("left", 0.0, 0.4), Source("left", 1.0, -0.1)], 1.75)
    assert [round(s, 6) for _, _, s in lines] == [0.4, 0.3]
    assert np.allclose(lines[1][0], [1.0, 1.5, 1.75]) and np.allclose(lines[1][1], [0.0, 2.0, 1.0])

def test_front_segments_right_corner_cancelled_by_control_jump():
    """Фронт из правого угла приходит к левому концу в момент t0 и гасится скачком управления той же величины."""
    lines = front_segments(P, [Source("right", 0.0, -0.35), Source("left", 0.5, -0.35)], 3.0)
    assert len(lines) == 1
    t, x, size = lines[0]
    assert np.allclose(t, [0.0, 0.5]) and np.allclose(x, [2.0, 0.0])
    lines = front_segments(P, [Source("left", 0.2, 0.3), Source("left", 0.7, 0.3)], 1.0)       # разные линии
    assert len(lines) == 2
