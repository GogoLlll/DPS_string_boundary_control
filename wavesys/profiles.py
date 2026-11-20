import numpy as np

from .segments import split_time

class ProfileSolver:
    def __init__(self, params, data, control):
        self.params, self.data, self.control = params, data, control

    def A0(self, s):
        x = self.params.a * np.asarray(s, dtype=float)
        return self.data.y0(x) + self.data.F(x) / self.params.a

    def B0(self, s):
        x = self.params.l - self.params.a * np.asarray(s, dtype=float)
        return self.data.y0(x) - self.data.F(x) / self.params.a

    def control_sum(self, k, s, parity):
        """Сумма u_j(s) по j < k с заданной чётностью (0 — чётные, 1 — нечётные)."""
        total = np.zeros(np.shape(s))
        for j in range(parity, k, 2):
            total = total + self.control.slot(j, s, self.params.t0)
        return total

    def segment_A(self, k, s):
        if k % 2 == 0:
            return self.A0(s) - 2 * self.control_sum(k, s, 0)
        return -self.B0(s) - 2 * self.control_sum(k, s, 1)

    def segment_B(self, k, s):
        if k % 2 == 0:
            return self.B0(s) + 2 * self.control_sum(k, s, 1)
        return -self.A0(s) + 2 * self.control_sum(k, s, 0)

    def _evaluate(self, z, segment):
        k, s = split_time(z, self.params.t0)
        out = np.zeros(np.shape(s))
        for kk in np.unique(k):
            mask = k == kk
            out[mask] = segment(int(kk), s[mask])
        return out

    def A(self, z):
        return self._evaluate(z, self.segment_A)

    def B(self, z):
        return self._evaluate(z, self.segment_B)
