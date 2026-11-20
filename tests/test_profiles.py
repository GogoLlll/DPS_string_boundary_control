import numpy as np

from wavesys.profiles import ProfileSolver

def test_initial_segments_reproduce_initial_data(smooth, inner_points):
    p, d = smooth.params, smooth.data
    pr = ProfileSolver(p, d, smooth.control)
    x = p.l * inner_points
    assert np.allclose(0.5 * (pr.A0(x / p.a) + pr.B0((p.l - x) / p.a)), d.y0(x), atol=1e-14)
    assert np.allclose(0.5 * p.a * (pr.A0(x / p.a) - pr.B0((p.l - x) / p.a)), d.F(x), atol=1e-14)

def test_segments_satisfy_recursion(smooth, inner_points):
    p = smooth.params
    pr = ProfileSolver(p, smooth.data, smooth.control)
    s = p.t0 * inner_points
    for k in range(9):
        assert np.allclose(pr.segment_A(k + 1, s), -pr.segment_B(k, s), atol=1e-13)
        assert np.allclose(pr.segment_B(k + 1, s), 2 * smooth.control.slot(k, s, p.t0) - pr.segment_A(k, s), atol=1e-13)

def test_first_steps_match_lecture_table(smooth, inner_points):
    """A2 = A0 − 2u0, B2 = B0 + 2u1, A3 = −B0 − 2u1, B3 = −A0 + 2u0 + 2u2 (п. 2.3.1)."""
    p = smooth.params
    pr = ProfileSolver(p, smooth.data, smooth.control)
    s = p.t0 * inner_points
    u = [smooth.control.slot(j, s, p.t0) for j in range(3)]
    assert np.allclose(pr.segment_A(2, s), pr.A0(s) - 2 * u[0])
    assert np.allclose(pr.segment_B(2, s), pr.B0(s) + 2 * u[1])
    assert np.allclose(pr.segment_A(3, s), -pr.B0(s) - 2 * u[1])
    assert np.allclose(pr.segment_B(3, s), -pr.A0(s) + 2 * u[0] + 2 * u[2])

def test_boundary_links(smooth):
    """A(t + t0) = −B(t), B(t + t0) = 2u(t) − A(t) — формула (2.10)."""
    p = smooth.params
    pr = ProfileSolver(p, smooth.data, smooth.control)
    t = np.linspace(0.01, 4.9, 333)
    assert np.allclose(pr.A(t + p.t0), -pr.B(t), atol=1e-13)
    assert np.allclose(pr.B(t + p.t0), 2 * smooth.control(t, p.t0) - pr.A(t), atol=1e-13)

def test_profile_keeps_array_shape(smooth):
    pr = ProfileSolver(smooth.params, smooth.data, smooth.control)
    assert pr.A(np.zeros((4, 3)) + 0.3).shape == (4, 3)
