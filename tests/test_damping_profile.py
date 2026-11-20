import numpy as np
import pytest

from wavesys.damping import BaseProfile
from wavesys.problem import Problem

Y0 = "0.6*(1 - x/l) + 0.7*sin(pi*x/l) + 0.2*sin(3*pi*x/l)"
Y1 = "0.5*sin(2*pi*x/l) + 0.3*sin(pi*x/l)"

@pytest.fixture(scope="module")
def generic():
    pr = Problem.from_text(l=1.5, a=2.0, T=1.0, y0=Y0, y1=Y1, u="0")
    return pr.params, pr.data

def lecture(kappa=0.8, l=1.0, a=1.0):
    pr = Problem.from_text(l=l, a=a, T=1.0, y0=f"{kappa}*(1 - x/l)*l", y1="0", u="0")     # y0 = kappa (l - x)
    return pr.params, pr.data

def test_base_profile_branches(generic):
    p, d = generic
    bp = BaseProfile(p, d, 2 * p.t0)
    s = np.linspace(0.01, p.t0 - 0.01, 17)
    assert np.allclose(bp.g(0, s), 0.5 * (d.y0(p.a * s) + d.F(p.a * s) / p.a))
    sigma = p.t0 + s
    assert np.allclose(bp.g(1, s), 0.5 * (-d.y0(2 * p.l - p.a * sigma) + d.F(2 * p.l - p.a * sigma) / p.a))

@pytest.mark.parametrize("q, expected", [
    (0.4, [(0.2, 1), (0.7, 0), (1.3, 0)]),          # T < t0
    (1.5, [(0.2, 1), (0.9, 1), (1.3, 1), (1.7, 0)]),
    (2.0, [(0.2, 1), (1.7, 1)]),
    (3.3, [(0.2, 2), (1.2, 2), (1.4, 1), (1.9, 1)]),  # rho = 1.3 t0
    (4.0, [(0.5, 2), (1.5, 2)]),
    (6.6, [(0.5, 4), (0.7, 3), (1.9, 3)]),
])
def test_repeat_count_formula_XX(generic, q, expected):
    p, d = generic
    bp = BaseProfile(p, d, q * p.t0)
    for sigma_frac, n in expected:
        half = int(sigma_frac >= 1)
        assert int(bp.count(half, np.array([(sigma_frac - half) * p.t0]))[0]) == n
    # прямой подсчёт моментов sigma + 2 m t0 < T
    for sigma_frac in (0.13, 0.77, 1.31, 1.93):
        direct = sum(1 for m in range(50) if sigma_frac + 2 * m < q - 1e-12)
        half = int(sigma_frac >= 1)
        assert int(bp.count(half, np.array([(sigma_frac - half) * p.t0]))[0]) == direct

def test_count_is_exact_when_T_is_multiple_of_t0_in_floating_point():
    pr = Problem.from_text(l=0.4, a=5.0, T=1.0, y0="x", y1="0", u="0")          # t0 = 0.08 не представимо точно
    for k in range(1, 40):
        bp = BaseProfile(pr.params, pr.data, k * pr.params.t0)
        assert (bp.kT, bp.sT) == (k, 0.0)
        assert int(bp.count(0, np.array([1e-9]))[0]) == (k + 1) // 2

def test_critical_time_reproduces_lecture_6(generic):
    p, d = generic
    bp = BaseProfile(p, d, 2 * p.t0)
    x = np.linspace(0, p.l, 200001)
    F = d.F(x)
    Fbar = np.trapezoid(F, x) / p.l
    assert bp.r_optimal() == pytest.approx(-Fbar / p.a, abs=1e-9)
    assert bp.cost_min() == pytest.approx(np.trapezoid(d.y0(x) ** 2 + (F - Fbar) ** 2 / p.a ** 2, x) / (2 * p.a), rel=1e-8)
    assert bp.cost_min() == pytest.approx(0.25582156, abs=5e-9)                 # reshenie.pdf, п. 3.3
    assert bp.r_optimal() == pytest.approx(-0.131303, abs=5e-7)

@pytest.mark.parametrize("q, cost, r", [(2.5, 0.21216354, -0.003806), (3.0, 0.17328330, 0.131389), (3.3, 0.16996482, 0.104710),
                                        (4.0, 0.12791078, -0.131303), (5.7, 0.09495204, -0.058343), (8.0, 0.06395539, -0.131303)])
def test_supercritical_numbers_from_the_written_solution(generic, q, cost, r):
    p, d = generic
    bp = BaseProfile(p, d, q * p.t0)
    assert bp.cost_min() == pytest.approx(cost, abs=5e-9)
    assert bp.r_optimal() == pytest.approx(r, abs=5e-7)

def test_cost_is_parabola_in_r_with_minimum_at_r_optimal(generic):
    p, d = generic
    bp = BaseProfile(p, d, 3.3 * p.t0)
    rs, e = bp.r_optimal(), 1e-3
    W = bp._weights()[0]
    assert bp.cost(rs) == pytest.approx(bp.cost_min(), abs=1e-14)
    assert bp.cost(rs + e) - bp.cost(rs) == pytest.approx(W / 4 * e * e, rel=1e-6)
    assert bp.cost(rs - e) - bp.cost(rs) == pytest.approx(W / 4 * e * e, rel=1e-6)

def test_multiples_of_two_passes_divide_the_cost(generic):
    p, d = generic
    base = BaseProfile(p, d, 2 * p.t0)
    for K in (2, 3, 5):
        bp = BaseProfile(p, d, 2 * K * p.t0)
        assert bp.cost_min() == pytest.approx(base.cost_min() / K, rel=1e-12)
        assert bp.r_optimal() == pytest.approx(base.r_optimal(), abs=1e-12)

def test_cost_does_not_increase_with_time(generic):
    p, d = generic
    costs = [BaseProfile(p, d, q * p.t0).cost_min() for q in np.arange(2.0, 8.01, 0.25)]
    assert all(b <= a + 1e-12 for a, b in zip(costs, costs[1:]))

@pytest.mark.parametrize("l, a", [(1.0, 1.0), (2.3, 0.7), (0.4, 5.0)])
def test_lecture_example_closed_forms(l, a):
    """y0 = κ(l − x): J* = κ²l³/(6a); J_min(3t0) = 11κ²l³/(96a), r*(3t0) = κl/6; J_min(4t0) = κ²l³/(12a)."""
    kappa = 0.8
    p, d = lecture(kappa, l, a)
    c = kappa ** 2 * l ** 3 / a
    assert BaseProfile(p, d, 2 * p.t0).cost_min() == pytest.approx(c / 6, rel=1e-10)
    assert BaseProfile(p, d, 2 * p.t0).r_optimal() == pytest.approx(0.0, abs=1e-12)
    three = BaseProfile(p, d, 3 * p.t0)
    assert three.cost_min() == pytest.approx(11 * c / 96, rel=1e-10)
    assert three.r_optimal() == pytest.approx(kappa * l / 6, rel=1e-10)
    assert BaseProfile(p, d, 4 * p.t0).cost_min() == pytest.approx(c / 12, rel=1e-10)

@pytest.mark.parametrize("q", [0.3, 0.5, 1.0, 1.5, 1.9])
def test_subcritical_residual_of_lecture_example(q):
    """Φ_min(T) = κ²a³(2t0 − T)³/24, r = −2·(среднее G) → y0(0) при T → 2t0."""
    kappa, l, a = 0.8, 1.5, 2.0
    p, d = lecture(kappa, l, a)
    bp = BaseProfile(p, d, q * p.t0)
    T = q * p.t0
    assert not bp.supercritical and bp.cost_min() is None
    assert bp.residual_min() == pytest.approx(kappa ** 2 * a ** 3 * (2 * p.t0 - T) ** 3 / 24, rel=1e-10)
    assert bp.r_optimal() == pytest.approx(kappa * a * T / 2, rel=1e-10)       # −2·(−κaT/4)

def test_subcritical_numbers_from_the_written_solution(generic):
    p, d = generic
    assert BaseProfile(p, d, 1.0 * p.t0).residual_min() == pytest.approx(0.05244892, abs=5e-9)
    assert BaseProfile(p, d, 1.7 * p.t0).residual_min() == pytest.approx(0.00337990, abs=5e-9)
    assert BaseProfile(p, d, 3.0 * p.t0).residual_min() == 0.0
    assert BaseProfile(p, d, 1.0 * p.t0).initial_norm() == pytest.approx(1.02328624, abs=5e-9)

def test_initial_norm_is_twice_a_times_critical_cost(generic):
    """Φ(0) = 2a·J* — связь из п. 3.6 решения."""
    p, d = generic
    bp = BaseProfile(p, d, 2 * p.t0)
    assert bp.initial_norm() == pytest.approx(2 * p.a * bp.cost_min(), rel=1e-12)

def test_minimal_time():
    generic = Problem.from_text(l=1.0, a=1.0, T=1.0, y0="0.8*(1 - x/l)", y1="0", u="0")
    assert BaseProfile(generic.params, generic.data, 1.0).minimal_time() == pytest.approx(2.0)
    rest = Problem.from_text(l=1.0, a=1.0, T=1.0, y0="0", y1="0", u="0")
    assert BaseProfile(rest.params, rest.data, 1.0).minimal_time() == 0.0
    wave = Problem.from_text(l=1.0, a=1.0, T=1.0, y0="0.6*cos(pi*x/1.2)^2*ind(x, 0, 0.6)",
                             y1="-0.6*pi/1.2*sin(pi*x/0.6)*ind(x, 0, 0.6)", u="0")           # y1 = a y0', носитель (0, 0.6)
    bp = BaseProfile(wave.params, wave.data, 0.6)
    assert bp.minimal_time() == pytest.approx(0.6, abs=2e-3)
    assert bp.residual_min() < 1e-20 and bp.r_optimal() == pytest.approx(0.6, abs=1e-10)
    two_waves = Problem.from_text(l=1.0, a=1.0, T=1.0, y0="0.5*sin(pi*x) + 0.7*sin(pi*(x - 0.45)/0.55)^2*ind(x, 0.45, 1)",
                                  y1="0.5*pi*cos(pi*x) - 0.7*pi/0.55*sin(2*pi*(x - 0.45)/0.55)*ind(x, 0.45, 1)", u="0")
    # волна влево на всём отрезке + волна вправо на (0.45, 1): T_min = (2l − d)/a = 1.55
    assert BaseProfile(two_waves.params, two_waves.data, 1.0).minimal_time() == pytest.approx(1.55, abs=2e-3)
