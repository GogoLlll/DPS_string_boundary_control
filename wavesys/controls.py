from abc import ABC, abstractmethod

import numpy as np

from .segments import before, split_time

class BoundaryControl(ABC):

    label = "управление"

    @abstractmethod
    def slot(self, j, s, t0):
        """Значения u(j t0 + s) для массива смещений s, 0 <= s < t0"""

    def __call__(self, t, t0):
        """Значения u(t) для массива моментов t >= 0"""
        k, s = split_time(t, t0)
        out = np.zeros(np.shape(s))
        for j in np.unique(k):
            mask = k == j
            out[mask] = self.slot(int(j), s[mask], t0)
        return out

    def jump_times(self, t0, horizon):
        return []

class ZeroControl(BoundaryControl):
    label = "u = 0"

    def slot(self, j, s, t0):
        return np.zeros(np.shape(s))

class FunctionControl(BoundaryControl):
    def __init__(self, f, active_until=None, label=None):
        self.f = f
        self.active_until = active_until
        self.label = label or getattr(f, "text", "u(t)")

    def slot(self, j, s, t0):
        s = np.asarray(s, dtype=float)
        values = np.asarray(self.f(j * t0 + s), dtype=float) + np.zeros(s.shape)
        if self.active_until is not None:
            k_end, s_end = split_time(self.active_until, t0)
            values = np.where(before(j, s, int(k_end), float(s_end)), values, 0.0)
        return values

    def jump_times(self, t0, horizon):
        times = find_jumps(lambda t: self(t, t0), 0.0, horizon)
        if self.active_until is not None and 0 < self.active_until < horizon:
            if not any(abs(t - self.active_until) < 1e-9 * max(horizon, 1.0) for t in times):
                left = float(np.asarray(self.f(np.array([self.active_until]))).ravel()[0])
                if abs(left) > 1e-12:
                    times.append(float(self.active_until))
        return sorted(times)

def find_jumps(f, lo, hi, samples=4001, rounds=6, fan=16):
    if not hi > lo:
        return []
    t = np.linspace(lo, hi, samples)
    v = np.asarray(f(t), dtype=float)
    scale = max(float(np.max(np.abs(v))), 1e-300)
    dv = np.abs(np.diff(v))
    cand = np.nonzero(dv > 1e-9 * scale)[0]
    a_, b_, d = t[cand], t[cand + 1], dv[cand]
    frac = np.linspace(0.0, 1.0, fan + 1)
    for _ in range(rounds):
        if a_.size == 0:
            break
        pts = a_[:, None] + (b_ - a_)[:, None] * frac[None, :]
        vals = np.asarray(f(pts.ravel()), dtype=float).reshape(pts.shape)
        inc = np.abs(np.diff(vals, axis=1))
        best = np.argmax(inc, axis=1)
        top = inc[np.arange(a_.size), best]
        keep = top > 0.8 * d
        rows = np.nonzero(keep)[0]
        a_, b_, d = pts[rows, best[rows]], pts[rows, best[rows] + 1], top[rows]
    jumps = []
    for b in np.sort(b_):
        if lo < b < hi and (not jumps or abs(b - jumps[-1]) > 1e-9 * (hi - lo)):
            jumps.append(float(b))
    return jumps
