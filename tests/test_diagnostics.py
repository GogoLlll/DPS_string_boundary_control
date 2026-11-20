import numpy as np

from wavesys import diagnostics
from wavesys.problem import Problem

def _by_name(checks):
    return {c.name: c for c in checks}

def test_all_checks_pass_for_smooth_problem(smooth):
    checks = diagnostics.run_all(smooth, smooth.solve(), fd_nodes=80)
    assert [c.name for c in checks if c.ok is False] == []
    assert sum(c.ok is True for c in checks) >= 6

def test_smoothness_verdicts():
    cases = {
        "классическое": dict(y0="0.5*(1 - x/l)", y1="0", u="0.5"),
        "вторые производные разрывны": dict(y0="0", y1="0", u="0.4*t^2"),
        "с изломами": dict(y0="0", y1="0", u="0.4*t"),
        "есть скачки": dict(y0="0.5*(1 - x/l)", y1="0", u="0"),
    }
    for expected, kw in cases.items():
        pr = Problem.from_text(l=1.0, a=1.0, T=2.0, **kw)
        verdict = _by_name(diagnostics.compatibility(pr, pr.solve()))["вывод о гладкости"].note
        assert expected in verdict, (expected, verdict)

def test_control_jump_makes_solution_generalized():
    pr = Problem.from_text(l=1.0, a=1.0, T=2.0, y0="0", y1="0", u="0.5*ind(t, 0.3, 0.6)")
    verdict = _by_name(diagnostics.compatibility(pr, pr.solve()))["вывод о гладкости"].note
    assert "есть скачки" in verdict

def test_checks_catch_a_wrong_solution(smooth):
    sol = smooth.solve()

    class Broken:
        profiles = sol.profiles
        state_norm = staticmethod(sol.state_norm)
        Q = staticmethod(sol.Q)
        grid = staticmethod(lambda t, x: sol.grid(t, x) * 1.01)
        y = staticmethod(lambda t, x: sol.y(t, x) + 0.01 * np.sin(7 * np.asarray(x) * np.asarray(t)))
        left_boundary = staticmethod(lambda t: sol.left_boundary(t) + 0.01)
        right_boundary = staticmethod(lambda t: sol.right_boundary(t) + 0.02)

    failed = {c.name for c in diagnostics.run_all(smooth, Broken, fd_nodes=40) if c.ok is False}
    assert "max |y(t,0) − u(t)|" in failed
    assert "max |y(t,l)|" in failed
    assert any("ромба" in n for n in failed)
    assert any("расхождение" in n for n in failed)

def test_fd_comparison_reports_front_nodes_for_discontinuous_solution():
    pr = Problem.from_text(l=1.0, a=1.0, T=2.0, y0="0.5*(1 - x/l)", y1="0", u="0")
    checks = diagnostics.fd_comparison(pr, pr.solve(), nodes=50)
    assert checks[0].ok is True
    assert len(checks) == 2 and checks[1].ok is None
