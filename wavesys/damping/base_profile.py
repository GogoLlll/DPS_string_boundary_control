import numpy as np

from ..segments import split_time

_GX, _GW = np.polynomial.legendre.leggauss(8)

class BaseProfile:
    def __init__(self, params, data, T, cells=256, breaks=()):
        self.params, self.data, self.T = params, data, float(T)
        k, s = split_time(self.T, params.t0)
        self.kT, self.sT = int(k), float(s)                 # T = kT t0 + sT
        self.K = self.kT // 2                               # T = 2 K t0 + rho
        self.rho = (self.kT % 2) * params.t0 + self.sT
        self.supercritical = self.kT >= 2                   # T >= 2 t0: в каждой сумме есть слагаемое
        t0 = params.t0
        marks = sorted({0.0, t0} | {float(b) for b in list(breaks) + [self.sT] if 0.0 < b < t0})
        nodes, weights = [], []
        for lo, hi in zip(marks[:-1], marks[1:]):
            n = max(1, int(round(cells * (hi - lo) / t0)))
            edges = np.linspace(lo, hi, n + 1)
            half = np.diff(edges) / 2
            mid = edges[:-1] + half
            nodes.append((mid[:, None] + half[:, None] * _GX[None, :]).ravel())
            weights.append((half[:, None] * _GW[None, :]).ravel())
        self.s, self.w = np.concatenate(nodes), np.concatenate(weights)

    def g(self, p, s):
        """G в точке sigma = p t0 + s."""
        pr, d = self.params, self.data
        s = np.asarray(s, dtype=float)
        if p == 0:
            x = pr.a * s
            return 0.5 * (d.y0(x) + d.F(x) / pr.a)
        x = pr.l - pr.a * s                                 # 2l - a sigma
        return 0.5 * (-d.y0(x) + d.F(x) / pr.a)

    def count(self, p, s):
        s = np.asarray(s, dtype=float)
        whole = max(0, (self.kT - p + 1) // 2)              # участки j = p, p+2, ... с j < kT
        partial = (s < self.sT) if self.kT % 2 == p else np.zeros(s.shape, dtype=bool)
        return whole + partial.astype(int)

    def slots(self, p):
        return list(range(p, self.kT + 1, 2))

    def integrate(self, f):
        """Интеграл по базовому интервалу (0, 2 t0) от f(p, s)."""
        return float(sum(np.sum(self.w * f(p, self.s)) for p in (0, 1)))

    def _weights(self):
        W = M = S = 0.0
        for p in (0, 1):
            n = self.count(p, self.s)
            on = n > 0
            inv = np.where(on, 1.0 / np.maximum(n, 1), 0.0)
            G = self.g(p, self.s)
            W += float(np.sum(self.w * inv))
            M += float(np.sum(self.w * inv * G))
            S += float(np.sum(self.w * inv * G * G))
        return W, M, S

    def tail(self):
        length = total = 0.0
        parts = []
        for p in (0, 1):
            off = self.count(p, self.s) == 0
            G = self.g(p, self.s)
            length += float(np.sum(self.w * off))
            total += float(np.sum(self.w * off * G))
            parts.append((off, G))
        if length <= 0:
            return 0.0, 0.0, 0.0
        mean = total / length
        spread = float(sum(np.sum(self.w * off * (G - mean) ** 2) for off, G in parts))
        return length, mean, spread

    def cost(self, r):
        W, M, S = self._weights()
        return S + r * M + r * r * W / 4

    def r_optimal(self):
        if self.supercritical:
            W, M, _ = self._weights()
            return -2 * M / W
        return -2 * self.tail()[1]

    def cost_min(self):
        if not self.supercritical:
            return None
        W, M, S = self._weights()
        return S - M * M / W

    def residual_min(self):
        return 2 * self.params.a * self.tail()[2]

    def initial_norm(self):
        mean = self.integrate(self.g) / (2 * self.params.t0)
        return 2 * self.params.a * self.integrate(lambda p, s: (self.g(p, s) - mean) ** 2)

    def sample(self, points=600):
        t0 = self.params.t0
        s = t0 * (np.arange(points) + 0.5) / points
        sigma = np.concatenate([s, t0 + s])
        G = np.concatenate([self.g(0, s), self.g(1, s)])
        n = np.concatenate([self.count(0, s), self.count(1, s)])
        return sigma, G, n

    def minimal_time(self, points=4000, tol=1e-6):
        sigma, G, _ = self.sample(points)
        scale = max(float(np.max(np.abs(G))), 1e-300)
        hi = np.maximum.accumulate(G[::-1])[::-1]
        lo = np.minimum.accumulate(G[::-1])[::-1]
        flat = (hi - lo) <= tol * scale
        if flat[0]:
            return 0.0
        last = int(np.nonzero(~flat)[0][-1])
        if last >= len(sigma) - 2:
            return 2 * self.params.t0
        step = sigma[1] - sigma[0]
        return float(min(sigma[last] + step, 2 * self.params.t0))
