import numpy as np

from ..controls import BoundaryControl, find_jumps
from ..segments import before
from .base_profile import BaseProfile

class _DampingControl(BoundaryControl):

    def __init__(self, params, data, T, r, cells=256):
        self.profile = BaseProfile(params, data, T, cells=cells)
        self.params, self.T, self.r = params, float(T), float(r)

    def _active(self, j, s):
        return before(j, np.asarray(s, dtype=float), self.profile.kT, self.profile.sT)

    def jump_times(self, t0, horizon):
        pr = self.profile
        cand = {k * t0 for k in range(1, pr.kT + 2)} | {k * t0 + pr.sT for k in range(0, pr.kT + 1)} | {self.T}
        cand |= set(find_jumps(lambda t: self(t, t0), 0.0, min(self.T, horizon)))
        eps = 1e-9 * max(t0, horizon)
        scale = max(1.0, float(np.max(np.abs(self(np.linspace(0, min(self.T, horizon), 801), t0)))))
        out = []
        for tau in sorted(cand):
            if not (eps < tau < horizon - eps):
                continue
            size = float(self(np.array([tau + eps]), t0)[0] - self(np.array([tau - eps]), t0)[0])
            if abs(size) > 1e-7 * scale and (not out or tau - out[-1] > 10 * eps):
                out.append(float(tau))
        return out

class EqualShareControl(_DampingControl):

    label = "успокоение: равные доли"

    def slot(self, j, s, t0):
        s = np.asarray(s, dtype=float)
        p = j % 2
        n = self.profile.count(p, s)
        value = (self.profile.g(p, s) + self.r / 2) / np.maximum(n, 1)
        return np.where(self._active(j, s) & (n > 0), value, 0.0)

class ContinuousControl(_DampingControl):

    label = "успокоение: непрерывное управление"

    def __init__(self, params, data, T, r, eps, cells=256):
        super().__init__(params, data, T, r, cells=cells)
        pr, t0 = self.profile, params.t0
        self.joints = pr.K if pr.rho > 0 else pr.K - 1      # число внутренних стыков t = 2 m t0
        if self.joints < 1:
            raise ValueError("Непрерывное управление по формуле (XLII) строится только при T > 2 t0")
        if not eps > 0:
            raise ValueError("Ширина переходного слоя должна быть положительной")
        self.eps = float(min(eps, pr.rho, 2 * t0 - pr.rho) if pr.rho > 0 else min(eps, t0))
        self.y00 = float(np.asarray(data.y0(np.array([0.0]))).ravel()[0])
        self.g_start = 0.5 * (self.y00 + self.r)            # G(0+) + r/2
        self.g_end = 0.5 * (-self.y00 + self.r)             # G(2 t0 -) + r/2
        self.g_rho = float(pr.g(pr.kT % 2, np.array([pr.sT]))[0] + self.r / 2) if pr.rho > 0 else 0.0
        self.v_joint = (self.g_start - self.y00) / self.joints
        self.n_start = pr.K + 1 if pr.rho > 0 else pr.K     # число повторов у sigma = 0+

    def _v(self, m):
        pr = self.profile
        if m == 0:
            return self.y00
        if pr.rho == 0 and m == pr.K:
            return 0.0
        return self.v_joint

    def slot(self, j, s, t0):
        s = np.asarray(s, dtype=float)
        pr, e = self.profile, self.eps
        p, m = j % 2, j // 2
        sigma = p * t0 + s
        n = pr.count(p, s)
        value = (pr.g(p, s) + self.r / 2) / np.maximum(n, 1)
        lam = lambda q: np.maximum(1.0 - q, 0.0)
        value = value + (self._v(m) - self.g_start / self.n_start) * lam(sigma / e)
        if m <= pr.K - 1:
            value = value + (self._v(m + 1) - self.g_end / pr.K) * lam((2 * t0 - sigma) / e)
        if pr.rho > 0:
            w_m = self.g_rho / pr.K if m < pr.K else 0.0
            value = value + np.where(sigma < pr.rho, (w_m - self.g_rho / (pr.K + 1)) * lam((pr.rho - sigma) / e), 0.0)
        return np.where(self._active(j, s) & (n > 0), value, 0.0)
