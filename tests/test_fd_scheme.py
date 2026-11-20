import numpy as np

from wavesys.fd_scheme import leapfrog
from wavesys.problem import Problem

def test_grid_layout(smooth):
    g = leapfrog(smooth.params, smooth.data, smooth.control, 50, 2.0)
    p = smooth.params
    assert g.x[0] == 0 and g.x[-1] == p.l and len(g.x) == 51
    assert g.t[-1] >= 2.0 and np.allclose(np.diff(g.t), p.l / 50 / p.a)
    assert g.Y.shape == (len(g.t), 51)

def test_scheme_agrees_with_profiles_for_smooth_problem(smooth):
    g = leapfrog(smooth.params, smooth.data, smooth.control, 120, 4 * smooth.params.t0)
    exact = smooth.solve().grid(g.t, g.x)
    assert np.max(np.abs(exact - g.Y)) < 1e-12

def test_scheme_agrees_for_other_parameters():
    pr = Problem.from_text(l=0.4, a=5.0, T=0.5, y0="0.2*x*(l - x)/l^2", y1="sin(pi*x/l)", u="0")
    g = leapfrog(pr.params, pr.data, pr.control, 80, pr.horizon)
    assert np.max(np.abs(pr.solve().grid(g.t, g.x) - g.Y)) < 1e-12

def test_boundary_values_come_from_control(smooth):
    g = leapfrog(smooth.params, smooth.data, smooth.control, 40, 2.0)
    assert np.allclose(g.Y[1:, 0], smooth.control(g.t[1:], smooth.params.t0), atol=1e-13)
    assert np.all(g.Y[:, -1] == 0)
