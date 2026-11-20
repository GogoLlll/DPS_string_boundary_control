import csv
from dataclasses import dataclass, field

import numpy as np

from . import diagnostics
from .problem import Problem
from .solution import Solution

@dataclass
class Results:
    problem: Problem
    solution: Solution
    t: np.ndarray                 # моменты времени сетки
    x: np.ndarray                 # точки сетки
    Y: np.ndarray                 # y(t_i, x_j)
    u: np.ndarray                 # управление в моменты t
    left: np.ndarray              # y(t, 0)
    right: np.ndarray             # y(t, l)
    z: np.ndarray                 # аргумент профилей
    A: np.ndarray
    B: np.ndarray
    norm_t: np.ndarray
    norm: np.ndarray              # Phi(t)
    energy: object
    sources: list = field(default_factory=list)
    checks: list = field(default_factory=list)

def compute(problem: Problem, nx=201, nt=241, fd_nodes=200, norm_points=121) -> Results:
    p = problem.params
    sol = problem.solve()
    t = np.linspace(0.0, problem.horizon, int(nt))
    x = np.linspace(0.0, p.l, int(nx))
    z = np.linspace(0.0, problem.horizon + p.t0, 1601)
    norm_t = np.linspace(0.0, problem.horizon, int(norm_points))
    sources = problem.sources()
    energy = None
    if not sources:
        energy = np.array([sol.energy(v) for v in norm_t])
    return Results(
        problem=problem, solution=sol, t=t, x=x, Y=sol.grid(t, x),
        u=problem.control(t, p.t0), left=sol.left_boundary(t), right=sol.right_boundary(t),
        z=z, A=sol.profiles.A(z), B=sol.profiles.B(z),
        norm_t=norm_t, norm=np.array([sol.state_norm(v) for v in norm_t]), energy=energy,
        sources=sources, checks=diagnostics.run_all(problem, sol, fd_nodes),
    )

def value_table(results: Results, rows=13, cols=11):
    p = results.problem.params
    t = np.linspace(0.0, results.problem.horizon, rows)
    x = np.linspace(0.0, p.l, cols)
    return t, x, results.solution.grid(t, x)

def export_csv(path, results: Results, rows=None, cols=None):
    if rows is None and cols is None:
        t, x, Y = results.t, results.x, results.Y
    else:
        t, x, Y = value_table(results, rows or 13, cols or 11)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["t \\ x"] + [f"{v:.10g}" for v in x])
        for ti, row in zip(t, Y):
            w.writerow([f"{ti:.10g}"] + [f"{v:.10g}" for v in row])
