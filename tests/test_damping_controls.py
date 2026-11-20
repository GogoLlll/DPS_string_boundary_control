import numpy as np
import pytest

from wavesys.damping import BaseProfile, ContinuousControl, EqualShareControl, linking
from wavesys.fd_scheme import leapfrog
from wavesys.problem import Problem

Y0 = "0.6*(1 - x/l) + 0.7*sin(pi*x/l) + 0.2*sin(3*pi*x/l)"
Y1 = "0.5*sin(2*pi*x/l) + 0.3*sin(pi*x/l)"

def data(l=1.5, a=2.0, y0=Y0, y1=Y1):
    pr = Problem.from_text(l=l, a=a, T=1.0, y0=y0, y1=y1, u="0")
    return pr.params, pr.data

def optimal(p, d, T):
    return EqualShareControl(p, d, T, BaseProfile(p, d, T).r_optimal())

@pytest.mark.parametrize("q", [2.0, 2.5, 3.0, 3.3, 4.0, 5.7, 31 ** 0.5, 8.0])
def test_optimal_control_brings_string_to_rest_and_stays(q):
    p, d = data()
    T = q * p.t0
    u = optimal(p, d, T)
    sol = Problem(p, d, u, T + 2 * p.t0).solve()
    norm0 = sol.state_norm(0.0)
    assert norm0 > 0.5
    assert sol.state_norm(T) < 1e-24 * norm0 + 1e-28
    assert sol.state_norm(T + 1.37 * p.t0) < 1e-24 * norm0 + 1e-28
    x = p.l * (np.arange(100) + 0.37) / 100
    assert np.max(np.abs(sol.y(T, x))) < 1e-13
    assert linking.cost(u, p, d, T) == pytest.approx(BaseProfile(p, d, T).cost_min(), rel=1e-9)

def test_state_is_not_at_rest_just_before_T():
    p, d = data()
    T = 3.3 * p.t0
    sol = Problem(p, d, optimal(p, d, T), T).solve()
    assert sol.state_norm(T - 0.2 * p.t0) > 1e-4

@pytest.mark.parametrize("l, a", [(1.0, 1.0), (2.3, 0.7), (0.4, 5.0), (7.0, 3.1)])
def test_rest_on_regular_grid_points_for_other_parameters(l, a):
    p, d = data(l=l, a=a)
    T = 3.3 * p.t0
    sol = Problem(p, d, optimal(p, d, T), T).solve()
    x = p.l * np.arange(1, 100) / 100
    assert np.max(np.abs(sol.y(T, x))) < 1e-12

def test_critical_control_is_lecture_6_formula():
    kappa, l, a = 0.8, 1.5, 2.0
    p, d = data(l, a, y0=f"{kappa}*(l - x)", y1="0")
    T = 2 * p.t0
    t = T * (np.arange(200) + 0.5) / 200
    h = optimal(p, d, T)
    assert np.allclose(h(t, p.t0), kappa / 2 * (l - a * t), atol=1e-13)
    u_en = EqualShareControl(p, d, T, kappa * l)
    assert np.allclose(u_en(t, p.t0), kappa * l - kappa * a * t / 2, atol=1e-13)
    assert linking.cost(u_en, p, d, T) == pytest.approx(4 * linking.cost(h, p, d, T), rel=1e-10)
    assert np.all(h(np.array([T, T + 0.1, T + 5.0]), p.t0) == 0.0)                # после T управление выключено

def test_three_passes_example_formula():
    kappa, l, a = 0.8, 1.0, 1.0
    p, d = data(l, a, y0=f"{kappa}*(l - x)", y1="0")
    u = optimal(p, d, 3 * p.t0)
    t1 = np.linspace(0.01, 0.99, 25)
    assert np.allclose(u(t1, p.t0), kappa / 4 * (l - a * t1) + kappa * l / 24, atol=1e-12)
    t2 = np.linspace(1.01, 1.99, 25)
    assert np.allclose(u(t2, p.t0), kappa / 2 * (l - a * t2) + kappa * l / 12, atol=1e-12)
    t3 = np.linspace(2.01, 2.99, 25)
    assert np.allclose(u(t3, p.t0), kappa / 4 * (l - a * (t3 - 2)) + kappa * l / 24, atol=1e-12)

def test_intermediate_states_shrink_linearly():
    p, d = data()
    K = 4
    T = 2 * K * p.t0
    sol = Problem(p, d, optimal(p, d, T), T).solve()
    x = p.l * (np.arange(60) + 0.37) / 60
    norm0 = sol.state_norm(0.0)
    for j in range(K + 1):
        assert np.allclose(sol.y(2 * j * p.t0, x), (1 - j / K) * d.y0(x), atol=1e-13)
        assert sol.state_norm(2 * j * p.t0) == pytest.approx((1 - j / K) ** 2 * norm0, rel=1e-9, abs=1e-25)

def test_any_r_gives_rest_but_costs_more():
    p, d = data()
    T = 3.3 * p.t0
    bp = BaseProfile(p, d, T)
    for r in (-0.8, 0.6, 2.0):
        u = EqualShareControl(p, d, T, r)
        assert Problem(p, d, u, T).solve().state_norm(T) < 1e-24
        assert linking.cost(u, p, d, T) == pytest.approx(bp.cost(r), rel=1e-9)
        assert bp.cost(r) > bp.cost_min()

@pytest.mark.parametrize("q", [0.3, 0.7, 1.0, 1.3, 1.7, 1.95])
def test_subcritical_control_reaches_minimal_residual(q):
    p, d = data()
    T = q * p.t0
    bp = BaseProfile(p, d, T)
    best = EqualShareControl(p, d, T, bp.r_optimal())
    reached = Problem(p, d, best, T).solve().state_norm(T)
    assert reached == pytest.approx(bp.residual_min(), rel=1e-7)
    assert reached > 0
    for r in (bp.r_optimal() + 0.5, bp.r_optimal() - 0.3):                       # другой r — остаток больше
        assert Problem(p, d, EqualShareControl(p, d, T, r), T).solve().state_norm(T) > reached * (1 + 1e-6)
    wiggle = Problem.from_text(l=1.5, a=2.0, T=T, y0=Y0, y1=Y1, u="0").solve().state_norm(T)   # без управления хуже
    assert wiggle > reached

def test_left_moving_wave_is_damped_before_critical_time():
    f = "0.9*cos(pi*x/(2*0.6*l))^2*ind(x, 0, 0.6*l)"
    df = "-a*0.9*pi/(2*0.6*l)*sin(pi*x/(0.6*l))*ind(x, 0, 0.6*l)"
    p, d = data(1.5, 2.0, y0=f, y1=df)
    T = 0.6 * p.t0
    u = optimal(p, d, T)
    t = T * (np.arange(100) + 0.5) / 100
    assert np.allclose(u(t, p.t0), d.y0(p.a * t), atol=1e-11)
    sol = Problem(p, d, u, T + p.t0).solve()
    assert sol.state_norm(0.0) > 0.05 and sol.state_norm(T) < 1e-20 and sol.state_norm(T + 0.5 * p.t0) < 1e-20
    x = p.l * (np.arange(50) + 0.37) / 50
    assert np.allclose(sol.y(0.3 * p.t0, x), d.y0(x + p.a * 0.3 * p.t0), atol=1e-11)
    assert Problem(p, d, optimal(p, d, 0.5 * p.t0), 0.5 * p.t0).solve().state_norm(0.5 * p.t0) > 1e-5   # раньше — нельзя

@pytest.mark.parametrize("q", [2.5, 3.3, 4.0, 5.7, 6.0])
@pytest.mark.parametrize("eps", [0.3, 0.05])
def test_continuous_control_is_continuous_exact_and_cheap(q, eps):
    p, d = data()
    T = q * p.t0
    bp = BaseProfile(p, d, T)
    u = ContinuousControl(p, d, T, bp.r_optimal(), eps * p.t0)
    t = np.linspace(0.0, T + 0.3 * p.t0, 60001)
    values = u(t, p.t0)
    assert values[0] == pytest.approx(float(d.y0(np.array([0.0]))[0]), abs=1e-12)        # u(0) = y0(0)
    slope_bound = 30 / (eps * p.t0)
    assert np.max(np.abs(np.diff(values))) < slope_bound * (t[1] - t[0])                 # нет скачков
    assert abs(float(u(np.array([T * (1 - 1e-9)]), p.t0)[0])) < 1e-6                     # u(T−) = 0
    prob = Problem(p, d, u, T + p.t0)
    assert prob.sources() == []
    sol = prob.solve()
    assert sol.state_norm(T) < 1e-24 and sol.state_norm(T + 0.7 * p.t0) < 1e-24
    cost = linking.cost(u, p, d, T)
    assert bp.cost_min() < cost < bp.cost_min() * (1 + 12 * eps)                         # близко к наименьшей

def test_continuous_cost_approaches_minimum_as_layer_shrinks():
    p, d = data()
    T = 4.0 * p.t0
    bp = BaseProfile(p, d, T)
    costs = [linking.cost(ContinuousControl(p, d, T, bp.r_optimal(), e * p.t0), p, d, T, cells=1024) for e in (0.3, 0.1, 0.03, 0.01)]
    assert all(b < a for a, b in zip(costs, costs[1:]))
    assert costs[-1] == pytest.approx(bp.cost_min(), rel=2e-2)
    assert costs[0] == pytest.approx(0.16789, abs=2e-5)                                  # reshenie.pdf, п. 3.8

def test_continuous_control_confirmed_by_independent_scheme():
    p, d = data()
    T = 3.3 * p.t0
    u = ContinuousControl(p, d, T, BaseProfile(p, d, T).r_optimal(), 0.1 * p.t0)
    g = leapfrog(p, d, u, 300, T + 2 * (p.l / 300) / p.a)
    n = int(round(T / (p.l / 300 / p.a)))
    inner = slice(1, -1)                                                                 # на концах значения задаёт граница
    assert np.max(np.abs(g.Y[n, inner])) < 1e-12                                         # покой в момент T
    assert np.max(np.abs(g.Y[n + 1, inner] - g.Y[n - 1, inner])) < 1e-12                 # и нулевая скорость

def test_continuous_control_needs_supercritical_time_and_positive_layer():
    p, d = data()
    with pytest.raises(ValueError):
        ContinuousControl(p, d, 2 * p.t0, 0.0, 0.1)
    with pytest.raises(ValueError):
        ContinuousControl(p, d, 1.5 * p.t0, 0.0, 0.1)
    with pytest.raises(ValueError):
        ContinuousControl(p, d, 3 * p.t0, 0.0, 0.0)

def test_layer_is_narrowed_to_fit():
    p, d = data()
    u = ContinuousControl(p, d, 2.1 * p.t0, 0.0, 5.0)                                    # rho = 0.1 t0
    assert u.eps == pytest.approx(0.1 * p.t0)

@pytest.mark.parametrize("q", [0.35, 0.8, 1.0, 1.6, 2.0, 3.3, 4.5])
def test_linking_residual_equals_state_norm_for_arbitrary_control(q):
    pr = Problem.from_text(l=1.5, a=2.0, T=q * 0.75, y0=Y0, y1=Y1, u="0.6*cos(3*t) + 0.4*sin(2*t)^2 - 0.3*sin(5*t + 1)")
    by_control, _ = linking.residual(pr.control, pr.params, pr.data, pr.horizon)
    assert by_control == pytest.approx(pr.solve().state_norm(pr.horizon), rel=1e-9)
    assert by_control > 1e-3

def test_linking_residual_without_control_is_initial_norm():
    pr = Problem.from_text(l=1.5, a=2.0, T=0.4, y0=Y0, y1=Y1, u="0")
    bp = BaseProfile(pr.params, pr.data, 1e-9)
    from wavesys.controls import ZeroControl
    assert linking.residual(ZeroControl(), pr.params, pr.data, 1e-9)[0] == pytest.approx(bp.initial_norm(), rel=1e-6)

def test_cost_of_simple_control():
    from wavesys.controls import FunctionControl
    p, d = data()
    u = FunctionControl(lambda t: 0.4 * t)
    T = 2.7 * p.t0
    assert linking.cost(u, p, d, T) == pytest.approx(0.16 * T ** 3 / 3, rel=1e-12)

def test_jump_times_of_optimal_control_cancel_in_pairs():
    from wavesys.characteristics import front_segments
    p, d = data()
    T = 3.3 * p.t0
    prob = Problem(p, d, optimal(p, d, T), T + p.t0)
    lines = front_segments(p, prob.sources(), prob.horizon)
    assert lines and all(times[-1] <= T + 1e-9 for times, _, _ in lines)                 # после T линий разрыва нет
