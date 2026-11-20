import io

import numpy as np

from .. import damping, design, results  # noqa: F401  (damping регистрирует построители и проверки)
from ..characteristics import front_segments, point_info
from ..functions import ALLOWED_FUNCTIONS, ExpressionError
from ..presets import DEFAULT, PRESETS
from ..problem import Problem
from .serialize import arr, num

class RequestError(ValueError):
    """Ошибка во входных данных"""

FORMULA_SOURCE = "Формула u(t)"

FUNCTION_HELP = [
    ("sin, cos, tan, exp, log, sqrt, abs", "обычные функции"),
    ("sinh, cosh, tanh, sign, min(a, b), max(a, b)", ""),
    ("pos(v)", "max(v, 0)"),
    ("ind(v, lo, hi)", "1 при lo ≤ v < hi, иначе 0"),
    ("step(v)", "1 при v ≥ 0, иначе 0"),
    ("where(условие, a, b)", "a там, где условие истинно, иначе b"),
    ("pi, e, l, a, t0, T", "постоянные и параметры задачи"),
    ("^ или **", "возведение в степень"),
]

def _number(payload, key, default=None, positive=False, integer=False, lo=None, hi=None):
    raw = payload.get(key, default)
    if raw is None or raw == "":
        raise RequestError(f"Не задано значение «{key}»")
    try:
        value = float(str(raw).replace(",", "."))
    except ValueError:
        raise RequestError(f"«{key}» должно быть числом, получено: {raw!r}") from None
    if not np.isfinite(value):
        raise RequestError(f"«{key}» должно быть конечным числом")
    if positive and not value > 0:
        raise RequestError(f"«{key}» должно быть положительным")
    if lo is not None and value < lo:
        raise RequestError(f"«{key}» не может быть меньше {lo:g}")
    if hi is not None and value > hi:
        raise RequestError(f"«{key}» не может быть больше {hi:g}")
    return int(round(value)) if integer else value

def meta():
    assert set(ALLOWED_FUNCTIONS) >= {"sin", "ind", "pos"}
    designers = [{"name": FORMULA_SOURCE, "description": "Управление задаётся формулой от t.", "options": {}}]
    for d in design.REGISTRY.values():
        designers.append({"name": d.name, "description": d.description,
                          "options": {k: {"label": v[0], "default": v[1], "choices": list(v[2]) if len(v) > 2 else None}
                                      for k, v in d.options.items()}})
    return {"presets": PRESETS, "default": DEFAULT, "sources": designers, "functions": FUNCTION_HELP}

def auto_frames(T, t0):
    return int(min(721, max(121, round(60 * T / t0) + 1)))

def build_problem(payload):
    l = _number(payload, "l", positive=True)
    a = _number(payload, "a", positive=True)
    T = _number(payload, "T", positive=True)
    y0, y1 = str(payload.get("y0", "0")), str(payload.get("y1", "0"))
    source = payload.get("source") or FORMULA_SOURCE
    if T / (l / a) > 400:
        raise RequestError("Слишком большое время: T не должно превышать 400 пробегов струны (T ≤ 400·l/a)")
    try:
        if source == FORMULA_SOURCE:
            return Problem.from_text(l, a, T, y0, y1, str(payload.get("u", "0")))
        if source not in design.REGISTRY:
            raise RequestError(f"Неизвестный источник управления: {source}")
        base = Problem.from_text(l, a, T, y0, y1, "0")
        options = dict(payload.get("options") or {})
        if not isinstance(payload.get("options") or {}, dict):
            raise RequestError("Параметры источника управления должны быть объектом")
        after = _number(options, "after", 0.0, lo=0.0, hi=50.0) if str(options.get("after", "")).strip() != "" else 0.0
        options.pop("after", None)
        control = design.REGISTRY[source].build(base.params, base.data, T, **options)
        texts = dict(base.texts, u=getattr(control, "label", source))
        return Problem(base.params, base.data, control, T + after * base.params.t0, texts=texts)
    except ExpressionError as exc:
        raise RequestError(str(exc)) from None
    except RequestError:
        raise
    except (ValueError, TypeError) as exc:
        raise RequestError(str(exc)) from None

def solve(payload):
    problem = build_problem(payload)
    p = problem.params
    nx = _number(payload, "nx", 201, integer=True, lo=11, hi=801)
    nt_raw = payload.get("nt")
    nt = auto_frames(problem.horizon, p.t0) if nt_raw in (None, "", "auto") else _number(payload, "nt", integer=True, lo=11, hi=2001)
    fd_nodes = _number(payload, "fd_nodes", 200, integer=True, lo=0, hi=1000)
    if fd_nodes and fd_nodes * fd_nodes * problem.horizon / p.t0 > 4e7:
        fd_nodes = max(20, int((4e7 * p.t0 / problem.horizon) ** 0.5))       # ограничение памяти и времени
    r = results.compute(problem, nx=nx, nt=nt, fd_nodes=fd_nodes, norm_points=min(nt, 241))
    return r, describe(r)

def describe(r):
    p, problem = r.problem.params, r.problem
    T = problem.horizon
    amp = float(np.nanmax(np.abs(r.Y))) if r.Y.size else 0.0
    passed = sum(c.ok is True for c in r.checks)
    judged = sum(c.ok is not None for c in r.checks)
    verdict = next((c.note for c in r.checks if c.name == "вывод о гладкости"), "")
    sources = [{"kind": s.kind, "time": num(s.time), "size": num(s.size), "label": s.label} for s in r.sources]
    fronts = [{"t": arr(pt), "x": arr(px), "size": num(size)} for pt, px, size in front_segments(p, r.sources, T)]
    return {
        "params": {"l": p.l, "a": p.a, "t0": p.t0, "T": T},
        "texts": problem.texts,
        "grid": {"t": arr(r.t), "x": arr(r.x), "Y": arr(r.Y)},
        "y0": arr(problem.data.y0(r.x)),
        "boundary": {"t": arr(r.t), "u": arr(r.u), "left": arr(r.left), "right": arr(r.right)},
        "profiles": {"z": arr(r.z), "A": arr(r.A), "B": arr(r.B)},
        "norm": {"t": arr(r.norm_t), "phi": arr(r.norm), "energy": None if r.energy is None else arr(r.energy)},
        "sources": sources,
        "fronts": fronts,
        "checks": [{"group": c.group, "name": c.name, "value": num(c.value), "limit": num(c.limit),
                    "ok": c.ok, "note": c.note} for c in r.checks],
        "damping": describe_damping(r),
        "stats": {
            "t0": p.t0, "T_over_t0": T / p.t0, "amplitude": amp,
            "phi_start": num(r.norm[0]), "phi_end": num(r.norm[-1]),
            "checks_passed": passed, "checks_total": judged, "verdict": verdict,
            "fronts": len(r.sources),
        },
    }

def describe_damping(r):
    """Отчёт о задаче успокоения — если управление построено модулем damping; иначе None."""
    info = getattr(r.problem.control, "damping", None)
    if info is None:
        return None
    T = info["target_time"]
    out = {k: (num(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else v) for k, v in info.items()
           if k not in ("sigma", "G", "count", "target")}
    out.update({"sigma": arr(info["sigma"]), "G": arr(info["G"]), "target": arr(info["target"]),
                "count": [int(v) for v in info["count"]],
                "state_norm_at_T": num(r.solution.state_norm(T)),
                "state_norm_at_end": num(r.solution.state_norm(r.problem.horizon))})
    return out

def point(r, payload):
    p = r.problem.params
    t = _number(payload, "t", lo=0.0)
    x = _number(payload, "x", lo=0.0, hi=p.l)
    info = point_info(p, t, x)
    pr = r.solution.profiles
    A = float(pr.segment_A(info["kA"], np.array([info["sA"]]))[0])
    B = float(pr.segment_B(info["kB"], np.array([info["sB"]]))[0])
    control = r.problem.control
    def moments(times):
        return [{"t": num(tau), "u": num(control(np.array([tau]), p.t0)[0])} for tau in times]
    return {
        "t": t, "x": x, "y": num(0.5 * (A + B)), "Q": num(0.5 * p.a * (A - B)),
        "A": num(A), "B": num(B),
        "zA": num(info["zA"]), "kA": info["kA"], "sA": num(info["sA"]),
        "zB": num(info["zB"]), "kB": info["kB"], "sB": num(info["sB"]),
        "A_from": info["A_from"], "B_from": info["B_from"], "region": info["region"],
        "control_A": moments(info["control_times_A"]), "control_B": moments(info["control_times_B"]),
    }

def _list(payload, key, lo, hi, limit=5):
    raw = payload.get(key, "")
    if isinstance(raw, (list, tuple)):
        raw = ",".join(str(v) for v in raw)
    values = []
    for piece in str(raw).replace(";", ",").split(","):
        piece = piece.strip()
        if not piece:
            continue
        try:
            v = float(piece)
        except ValueError:
            raise RequestError(f"«{piece}» не является числом") from None
        if not (lo <= v <= hi):
            raise RequestError(f"Значение {v:g} вне допустимого интервала [{lo:g}, {hi:g}]")
        values.append(v)
    if not values:
        raise RequestError("Не задано ни одного значения")
    if len(values) > limit:
        raise RequestError(f"Можно задать не больше {limit} значений")
    return sorted(set(values))


def sections(r, payload):
    p, T = r.problem.params, r.problem.horizon
    times = _list(payload, "times", 0.0, T)
    points = _list(payload, "points", 0.0, p.l)
    x = np.linspace(0.0, p.l, 601)
    t = np.linspace(0.0, T, 1201)
    sol = r.solution
    return {
        "x": arr(x), "shapes": [{"t": v, "y": arr(sol.y(np.full(x.shape, v), x))} for v in times],
        "t": arr(t), "traces": [{"x": v, "y": arr(sol.y(t, np.full(t.shape, v)))} for v in points],
    }

def table(r, payload):
    rows = _number(payload, "rows", 13, integer=True, lo=2, hi=201)
    cols = _number(payload, "cols", 11, integer=True, lo=2, hi=101)
    t, x, Y = results.value_table(r, rows, cols)
    return {"t": arr(t), "x": arr(x), "Y": arr(Y)}

def table_csv(r, payload):
    rows = _number(payload, "rows", 13, integer=True, lo=2, hi=2001)
    cols = _number(payload, "cols", 11, integer=True, lo=2, hi=801)
    import csv
    t, x, Y = results.value_table(r, rows, cols)
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(["t \\ x"] + [f"{v:.10g}" for v in x])
    for ti, row in zip(t, Y):
        w.writerow([f"{ti:.10g}"] + [f"{v:.10g}" for v in row])
    return buf.getvalue()
