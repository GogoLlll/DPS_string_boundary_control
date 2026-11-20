from dataclasses import dataclass

import numpy as np

from .characteristics import find_sources, front_position
from .fd_scheme import leapfrog

@dataclass
class Check:
    group: str
    name: str
    value: float
    limit: float          # допустимое значение; nan — строка справочная
    ok: object            # True / False / None (справочная строка)
    note: str = ""

def _scale(problem, solution):
    p = problem.params
    x = np.linspace(0, p.l, 201)
    t = np.linspace(0, problem.horizon, 201)
    return max(1.0, float(np.max(np.abs(problem.data.y0(x)))), float(np.max(np.abs(problem.control(t, p.t0)))),
               float(np.max(np.abs(problem.data.F(x)))) / p.a)

def initial_conditions(problem, solution):
    p = problem.params
    x = p.l * (np.arange(400) + 0.5) / 400
    tol = 1e-9 * _scale(problem, solution)
    dy = float(np.max(np.abs(solution.y(0.0, x) - problem.data.y0(x))))
    dq = solution.Q(0.0, x) - problem.data.F(x)                 # должно быть постоянной: Q_x = y_t = y1
    spread = float(np.max(dq) - np.min(dq))
    return [
        Check("Начальные условия", "max |y(0,x) − y0(x)|", dy, tol, dy <= tol),
        Check("Начальные условия", "размах Q(0,x) − F(x)", spread, tol, spread <= tol,
              "постоянство означает y_t(0,x) = y1(x)"),
    ]

def boundary_conditions(problem, solution):
    p = problem.params
    t = problem.horizon * (np.arange(600) + 0.5) / 600
    tol = 1e-9 * _scale(problem, solution)
    left = float(np.max(np.abs(solution.left_boundary(t) - problem.control(t, p.t0))))
    right = float(np.max(np.abs(solution.right_boundary(t))))
    return [
        Check("Граничные условия", "max |y(t,0) − u(t)|", left, tol, left <= tol),
        Check("Граничные условия", "max |y(t,l)|", right, tol, right <= tol),
    ]

def rhombus_identity(problem, solution, samples=3000, seed=0):
    p = problem.params
    rng = np.random.default_rng(seed)
    dt = rng.uniform(0.002, 0.2, samples) * min(p.t0, problem.horizon / 2)
    dx = p.a * dt
    x = rng.uniform(dx, p.l - dx)
    t = rng.uniform(dt, np.maximum(problem.horizon - dt, dt * 1.0001))
    r = solution.y(t + dt, x) + solution.y(t - dt, x) - solution.y(t, x + dx) - solution.y(t, x - dx)
    tol = 1e-9 * _scale(problem, solution)
    worst = float(np.max(np.abs(r)))
    return [Check("Уравнение", f"тождество ромба, {samples} случайных ромбов", worst, tol, worst <= tol,
                  "равносильно y_tt = a² y_xx, не требует гладкости")]

def compatibility(problem, solution):
    p, d = problem.params, problem.data
    h = 1e-4 * p.l
    ht = h / p.a
    f = lambda fun, pts: np.asarray(fun(np.array(pts, dtype=float)), dtype=float)
    y0 = f(d.y0, [0, h, 2 * h, 3 * h, p.l, p.l - h, p.l - 2 * h, p.l - 3 * h])
    y1 = f(d.y1, [0.0, p.l])
    u = problem.control(np.array([0, ht, 2 * ht, 3 * ht]), p.t0)
    du0 = (-3 * u[0] + 4 * u[1] - u[2]) / (2 * ht)
    d2u0 = (2 * u[0] - 5 * u[1] + 4 * u[2] - u[3]) / ht ** 2
    d2y0_left = (2 * y0[0] - 5 * y0[1] + 4 * y0[2] - y0[3]) / h ** 2
    d2y0_right = (2 * y0[4] - 5 * y0[5] + 4 * y0[6] - y0[7]) / h ** 2
    scale = _scale(problem, solution)
    rows = [
        ("порядок 0: u(0) − y0(0)", u[0] - y0[0], 1e-9 * scale),
        ("порядок 0: y0(l)", y0[4], 1e-9 * scale),
        ("порядок 1: u′(0) − y1(0)", du0 - y1[0], 1e-5 * scale / p.t0),
        ("порядок 1: y1(l)", y1[1], 1e-9 * scale / p.t0),
        ("порядок 2: u″(0) − a² y0″(0)", d2u0 - p.a ** 2 * d2y0_left, 1e-2 * scale / p.t0 ** 2),
        ("порядок 2: a² y0″(l)", p.a ** 2 * d2y0_right, 1e-2 * scale / p.t0 ** 2),
    ]
    checks = [Check("Согласование (2.2)–(2.4)", n, float(v), float(lim), None,
                    "выполнено" if abs(v) <= lim else "нарушено") for n, v, lim in rows]
    ok = [abs(v) <= lim for _, v, lim in rows]
    jumps = [s for s in problem.sources() if s.kind == "left" and s.time > 0]
    if not (ok[0] and ok[1]) or jumps:
        verdict = "решение обобщённое: есть скачки вдоль характеристик"
    elif not (ok[2] and ok[3]):
        verdict = "решение непрерывно, но с изломами от углов"
    elif not (ok[4] and ok[5]):
        verdict = "решение непрерывно дифференцируемо; вторые производные разрывны"
    else:
        verdict = "условия согласования выполнены: решение классическое (теорема 2.2)"
    checks.append(Check("Согласование (2.2)–(2.4)", "вывод о гладкости", float("nan"), float("nan"), None, verdict))
    return checks

def fd_comparison(problem, solution, nodes=200):
    p = problem.params
    grid = leapfrog(p, problem.data, problem.control, nodes, problem.horizon)
    exact = solution.grid(grid.t, grid.x)
    diff = np.abs(exact - grid.Y)
    dx = p.l / nodes
    on_front = np.zeros(diff.shape, dtype=bool)
    sources = find_sources(p, problem.data, problem.control, grid.t[-1])
    for src in sources:
        pos = front_position(p, src, grid.t)
        on_front |= np.abs(grid.x[None, :] - pos[:, None]) <= 0.51 * dx
    tol = 1e-9 * _scale(problem, solution)
    off = float(np.max(diff[~on_front])) if np.any(~on_front) else 0.0
    checks = [Check("Независимая схема", f"max расхождение в узлах вне фронтов ({nodes} отрезков)", off, tol, off <= tol,
                    "схема «крест» точна в узлах при a·dt = dx")]
    if sources:
        differing = int(np.count_nonzero(diff[on_front] > tol))
        total = int(np.count_nonzero(on_front))
        checks.append(Check("Независимая схема", "узлов на фронтах с другим значением", float(differing), float("nan"), None,
                            f"из {total} узлов на линиях разрыва; там значение зависит от выбранной стороны скачка"))
    return checks

EXTRA_CHECKS: list = []

def run_all(problem, solution, fd_nodes=200):
    checks = []
    checks += compatibility(problem, solution)
    checks += initial_conditions(problem, solution)
    checks += boundary_conditions(problem, solution)
    checks += rhombus_identity(problem, solution)
    if fd_nodes:
        checks += fd_comparison(problem, solution, fd_nodes)
    for extra in EXTRA_CHECKS:
        checks += extra(problem, solution)
    return checks
