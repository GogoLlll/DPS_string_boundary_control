from .. import diagnostics
from ..diagnostics import Check
from . import linking

def _breaks(problem):
    t0 = problem.params.t0
    return [s.time % t0 for s in problem.sources() if s.kind == "left" and s.time > 0]

def linking_consistency(problem, solution):
    p, H = problem.params, problem.horizon
    by_control, _ = linking.residual(problem.control, p, problem.data, H, breaks=_breaks(problem))
    by_state = solution.state_norm(H)
    scale = max(by_control, by_state, solution.state_norm(0.0), 1e-300)
    tol = (5e-3 if problem.sources() else 1e-7) * scale
    diff = abs(by_control - by_state)
    return [Check("Условие связи", "|Φ по управлению − Φ по состоянию| в конце наблюдения", diff, tol, diff <= tol,
                  "одно число получено только из управления, другое — только из движения струны")]

def damping_goal(problem, solution):
    info = getattr(problem.control, "damping", None)
    if info is None:
        return []
    T = info["target_time"]
    phi = solution.state_norm(T)
    norm0 = max(info["initial_norm"], 1e-300)
    group = "Успокоение"
    if info["exact"]:
        tol = 1e-9 * norm0
        rows = [Check(group, "норма состояния Φ(T) в момент успокоения", phi, tol, phi <= tol, "ноль означает покой")]
        if problem.horizon > T * (1 + 1e-12):
            later = solution.state_norm(problem.horizon)
            rows.append(Check(group, "норма состояния в конце наблюдения", later, tol, later <= tol,
                              "после T управление выключено, покой должен сохраняться"))
        return rows
    tol = 1e-6 * norm0
    diff = abs(phi - info["residual"])
    return [
        Check(group, "остаток Φ(T)", phi, float("nan"), None,
              "точное успокоение за это время невозможно; показан достигнутый остаток"),
        Check(group, "|Φ(T) − остаток, предсказанный по управлению|", diff, tol, diff <= tol),
        Check(group, "наименьший возможный остаток Φ_min(T)", info["residual_min"], float("nan"), None,
              "формула (XXXV); совпадает с достигнутым, если параметр r наилучший"),
    ]

for _check in (linking_consistency, damping_goal):
    if _check not in diagnostics.EXTRA_CHECKS:
        diagnostics.EXTRA_CHECKS.append(_check)
