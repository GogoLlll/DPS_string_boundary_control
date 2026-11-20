import numpy as np

from .profiles import ProfileSolver

_GAUSS_X, _GAUSS_W = np.polynomial.legendre.leggauss(8)

class Solution:
    def __init__(self, params, data, control):
        self.params, self.data, self.control = params, data, control
        self.profiles = ProfileSolver(params, data, control)

    def y(self, t, x):
        t, x = np.broadcast_arrays(np.asarray(t, dtype=float), np.asarray(x, dtype=float))
        p = self.params
        return 0.5 * (self.profiles.A(t + x / p.a) + self.profiles.B(t + (p.l - x) / p.a))

    def Q(self, t, x):
        t, x = np.broadcast_arrays(np.asarray(t, dtype=float), np.asarray(x, dtype=float))
        p = self.params
        return 0.5 * p.a * (self.profiles.A(t + x / p.a) - self.profiles.B(t + (p.l - x) / p.a))

    def grid(self, t, x):
        return self.y(np.asarray(t, dtype=float)[:, None], np.asarray(x, dtype=float)[None, :])

    def left_boundary(self, t):
        return self.y(t, np.zeros(np.shape(t)))

    def right_boundary(self, t):
        return self.y(t, np.full(np.shape(t), self.params.l))

    def _quadrature(self, t, pieces=96):
        t0 = self.params.t0
        joint = (np.floor(t / t0) + 1) * t0
        bounds = [t, joint, t + t0] if t < joint < t + t0 else [t, t + t0]
        zs, ws = [], []
        for lo, hi in zip(bounds[:-1], bounds[1:]):
            n = max(1, int(round(pieces * (hi - lo) / t0)))
            edges = np.linspace(lo, hi, n + 1)
            half = np.diff(edges) / 2
            mid = edges[:-1] + half
            zs.append((mid[:, None] + half[:, None] * _GAUSS_X[None, :]).ravel())
            ws.append((half[:, None] * _GAUSS_W[None, :]).ravel())
        return np.concatenate(zs), np.concatenate(ws)

    def state_norm(self, t):
        p = self.params
        z, w = self._quadrature(float(t))
        A, B = self.profiles.A(z), self.profiles.B(z)
        c = float(np.sum(w * (A - B) / 2) / p.t0)
        return float(p.a / 2 * np.sum(w * ((A - c) ** 2 + (B + c) ** 2)))

    def energy(self, t, h_rel=1e-6):
        p = self.params
        z, w = self._quadrature(float(t))
        h = h_rel * p.t0
        dA = (self.profiles.A(z + h) - self.profiles.A(z - h)) / (2 * h)
        dB = (self.profiles.B(z + h) - self.profiles.B(z - h)) / (2 * h)
        return float(np.sum(w * (dA ** 2 + dB ** 2)) / (4 * p.a))
