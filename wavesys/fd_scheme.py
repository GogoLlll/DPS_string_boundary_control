from dataclasses import dataclass

import numpy as np


@dataclass
class GridSolution:
    t: np.ndarray        # моменты слоёв, форма (nT + 1,)
    x: np.ndarray        # узлы, форма (N + 1,)
    Y: np.ndarray        # значения, форма (nT + 1, N + 1)

def leapfrog(params, data, control, nodes, horizon):
    N = int(nodes)
    l, a, t0 = params.l, params.a, params.t0
    dx = l / N
    dt = dx / a
    steps = max(1, int(np.ceil(horizon / dt - 1e-9)))
    x = dx * np.arange(N + 1)
    Y = np.zeros((steps + 1, N + 1))
    Y[0] = data.y0(x)
    Y[0, N] = 0.0
    n = np.arange(steps + 1)
    u = np.zeros(steps + 1)
    for j in np.unique(n // N):
        mask = (n // N) == j
        u[mask] = control.slot(int(j), (n[mask] % N) * dt, t0)
    F = data.F(x)
    Y[1, 1:N] = 0.5 * (Y[0, 2:] + Y[0, :-2]) + (F[2:] - F[:-2]) / (2 * a)
    Y[1, 0] = u[1]
    for m in range(1, steps):
        Y[m + 1, 1:N] = Y[m, 2:] + Y[m, :-2] - Y[m - 1, 1:N]
        Y[m + 1, 0] = u[m + 1]
    return GridSolution(t=dt * n, x=x, Y=Y)
