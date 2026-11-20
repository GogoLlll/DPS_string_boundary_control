import numpy as np
import pytest

from wavesys.problem import Problem

def test_initial_and_boundary_conditions(smooth, inner_points):
    p, sol = smooth.params, smooth.solve()
    x = p.l * inner_points
    assert np.allclose(sol.y(0.0, x), smooth.data.y0(x), atol=1e-14)
    t = 5 * p.t0 * inner_points
    assert np.allclose(sol.left_boundary(t), smooth.control(t, p.t0), atol=1e-13)
    assert np.allclose(sol.right_boundary(t), 0.0, atol=1e-13)

def test_initial_velocity(smooth, inner_points):
    p, sol = smooth.params, smooth.solve()
    x = p.l * inner_points[2:-2]
    h = 1e-6
    yt = (sol.y(h, x) - sol.y(0.0, x)) / h
    assert np.allclose(yt, smooth.data.y1(x), atol=1e-4)

def test_q_is_primitive_of_velocity(smooth, inner_points):
    p, sol = smooth.params, smooth.solve()
    x = p.l * inner_points[5:-5]
    t, h = 1.37 * p.t0, 1e-5
    qx = (sol.Q(t, x + h) - sol.Q(t, x - h)) / (2 * h)
    yt = (sol.y(t + h / p.a, x) - sol.y(t - h / p.a, x)) / (2 * h / p.a)
    assert np.allclose(qx, yt, atol=1e-5)

def test_standing_wave_example_1_1():
    pr = Problem.from_text(l=2.0, a=3.0, T=3.0, y0="0.5*sin(2*pi*x/l)", y1="0", u="0")
    sol, p = pr.solve(), pr.params
    t = np.linspace(0, 3, 41)[:, None]
    x = np.linspace(0, 2, 53)[None, :]
    exact = 0.5 * np.sin(2 * np.pi * x / p.l) * np.cos(2 * np.pi * p.a * t / p.l)
    assert np.allclose(sol.y(t, x), exact, atol=1e-13)

def test_cauchy_region_formula_2_18(smooth):
    p, d, sol = smooth.params, smooth.data, smooth.solve()
    t = 0.3 * p.t0
    x = np.linspace(p.a * t + 0.01, p.l - p.a * t - 0.01, 31)
    exact = 0.5 * (d.y0(x + p.a * t) + d.y0(x - p.a * t)) + (d.F(x + p.a * t) - d.F(x - p.a * t)) / (2 * p.a)
    assert np.allclose(sol.y(t, x), exact, atol=1e-14)

def test_first_layer_regions_II_III_IV(smooth):
    p, d, sol = smooth.params, smooth.data, smooth.solve()
    l, a, t0 = p.l, p.a, p.t0
    u = lambda t: smooth.control(np.asarray(t, dtype=float), t0)
    # II: x/a < t < (l − x)/a
    t, x = 0.5 * t0, 0.25 * l
    f2 = u(t - x / a) + 0.5 * (d.y0(x + a * t) - d.y0(a * t - x)) + (d.F(x + a * t) - d.F(a * t - x)) / (2 * a)
    assert sol.y(t, x) == pytest.approx(float(f2), abs=1e-13)
    # III: (l − x)/a < t < x/a
    t, x = 0.5 * t0, 0.75 * l
    f3 = 0.5 * (d.y0(x - a * t) - d.y0(2 * l - x - a * t)) + (d.F(2 * l - x - a * t) - d.F(x - a * t)) / (2 * a)
    assert sol.y(t, x) == pytest.approx(float(f3), abs=1e-13)
    # IV: t > max(x/a, (l − x)/a)
    t, x = 0.8 * t0, 0.45 * l
    f4 = (u(t - x / a) - 0.5 * (d.y0(2 * l - x - a * t) + d.y0(a * t - x))
          + (d.F(2 * l - x - a * t) - d.F(a * t - x)) / (2 * a))
    assert sol.y(t, x) == pytest.approx(float(f4), abs=1e-13)

def test_example_2_1_linear_end_motion():
    V, l, a = 0.4, 2.0, 0.5
    pr = Problem.from_text(l=l, a=a, T=2 * l / a, y0="0", y1="0", u=f"{V}*t")
    sol, t0 = pr.solve(), l / a
    assert sol.y(t0 / 2, l / 4) == pytest.approx(V * t0 / 4, abs=1e-13)
    assert sol.y(3 * t0 / 2, 3 * l / 4) == pytest.approx(V * t0 / 2, abs=1e-13)
    x = 0.3 * l
    for t in np.linspace(0.01, 1.99, 37) * t0:
        if t < x / a:
            exact = 0.0
        elif t < (2 * l - x) / a:
            exact = V * (t - x / a)
        else:
            exact = 2 * V * (l - x) / a
        assert sol.y(t, x) == pytest.approx(exact, abs=1e-13)

def test_response_series_formula_2_22():
    l, a = 1.0, 1.0
    pr = Problem.from_text(l=l, a=a, T=7.0, y0="0", y1="0", u="sin(3*t) + 0.2*t")
    sol, t0 = pr.solve(), l / a
    uf = lambda s: np.where(s > 0, np.sin(3 * s) + 0.2 * s, 0.0)
    t = np.linspace(0.05, 6.9, 61)[:, None]
    x = np.linspace(0.03, 0.97, 23)[None, :]
    q = sum(uf(t - x / a - 2 * m * t0) - uf(t - (2 * l - x) / a - 2 * m * t0) for m in range(6))
    assert np.allclose(sol.y(t, x), q, atol=1e-13)

def test_example_2_4_stationary_solution():
    pr = Problem.from_text(l=1.0, a=1.0, T=5.0, y0="0.5*(1 - x/l)", y1="0", u="0.5")
    sol = pr.solve()
    t = np.linspace(0, 5, 31)[:, None]
    x = np.linspace(0, 1, 27)[None, :]
    assert np.allclose(sol.y(t, x), 0.5 * (1 - x) + 0 * t, atol=1e-14)

def test_example_2_4_release_creates_jump():
    H = 0.5
    pr = Problem.from_text(l=1.0, a=1.0, T=1.0, y0=f"{H}*(1 - x/l)", y1="0", u="0")
    sol = pr.solve()
    t, x = 0.6, np.array([0.2, 0.5, 0.7, 0.9])
    exact = H * (1 - x) - H * (x < t)
    assert np.allclose(sol.y(t, x), exact, atol=1e-14)

def test_independent_of_primitive_constant(smooth, inner_points):
    from wavesys.initial_data import InitialData
    from wavesys.solution import Solution
    p, d = smooth.params, smooth.data
    shifted = InitialData(d.y0, d.y1, lambda x: d.F(x) + 0.83)
    a, b = smooth.solve(), Solution(p, shifted, smooth.control)
    t = 4 * p.t0 * inner_points[:, None]
    x = p.l * inner_points[None, ::3]
    assert np.allclose(a.y(t, x), b.y(t, x), atol=1e-13)

def test_grid_shape(smooth):
    assert smooth.solve().grid(np.linspace(0, 1, 7), np.linspace(0, 1.5, 5)).shape == (7, 5)

def test_state_norm_at_start_equals_norm_of_initial_data(smooth):
    p, d, sol = smooth.params, smooth.data, smooth.solve()
    x = np.linspace(0, p.l, 20001)
    F = d.F(x)
    integrand = d.y0(x) ** 2 + (F - np.trapezoid(F, x) / p.l) ** 2 / p.a ** 2
    assert sol.state_norm(0.0) == pytest.approx(np.trapezoid(integrand, x), rel=1e-7)

def test_state_norm_is_zero_at_rest_and_positive_otherwise():
    rest = Problem.from_text(l=1.0, a=1.0, T=1.0, y0="0", y1="0", u="0").solve()
    assert rest.state_norm(0.3) == 0.0
    moving = Problem.from_text(l=1.0, a=1.0, T=1.0, y0="sin(pi*x)", y1="0", u="0").solve()
    assert moving.state_norm(0.3) > 0.1

def test_state_norm_vanishes_under_lecture_6_control():
    l, a = 1.5, 2.0
    y0 = "0.6*(1 - x/l) + 0.7*sin(pi*x/l)"
    G1 = "0.5*((0.6*(1 - a*t/l) + 0.7*sin(pi*a*t/l)) + 0.6)"
    G2 = "0.5*(-(0.6*(1 - (2*l - a*t)/l) + 0.7*sin(pi*(2*l - a*t)/l)) + 0.6)"
    u = f"({G1})*ind(t, 0, t0) + ({G2})*ind(t, t0, 2*t0)"
    pr = Problem.from_text(l=l, a=a, T=3.0, y0=y0, y1="0", u=u)
    sol, t0 = pr.solve(), l / a
    assert sol.state_norm(0.0) > 0.1
    assert sol.state_norm(2 * t0) < 1e-24
    assert sol.state_norm(2.7 * t0) < 1e-24
    assert np.max(np.abs(sol.y(2 * t0, np.linspace(0.01, l - 0.01, 50)))) < 1e-13

def test_energy_is_conserved_without_control():
    pr = Problem.from_text(l=1.0, a=2.0, T=3.0, y0="0.5*sin(pi*x/l) + 0.1*sin(2*pi*x/l)", y1="0.3*sin(pi*x/l)", u="0")
    sol = pr.solve()
    e = [sol.energy(t) for t in np.linspace(0, 3, 13)]
    assert e[0] > 0.1
    assert np.allclose(e, e[0], rtol=1e-6)

def test_energy_of_standing_wave():
    H, l = 0.5, 2.0
    sol = Problem.from_text(l=l, a=3.0, T=1.0, y0=f"{H}*sin(pi*x/l)", y1="0", u="0").solve()
    assert sol.energy(0.37) == pytest.approx(H ** 2 * np.pi ** 2 / (4 * l), rel=1e-6)
