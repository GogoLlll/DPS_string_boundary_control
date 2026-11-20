import numpy as np

from ..segments import before
from .base_profile import BaseProfile

def _slot_values(control, profile, j):
    t0 = profile.params.t0
    values = np.asarray(control.slot(j, profile.s, t0), dtype=float)
    return np.where(before(j, profile.s, profile.kT, profile.sT), values, 0.0)

def stacked_sum(control, profile, p):
    total = np.zeros(profile.s.shape)
    for j in profile.slots(p):
        total = total + _slot_values(control, profile, j)
    return total

def residual(control, params, data, T, cells=256, breaks=()):
    profile = BaseProfile(params, data, T, cells=cells, breaks=breaks)
    D = [stacked_sum(control, profile, p) - profile.g(p, profile.s) for p in (0, 1)]
    mean = float(sum(np.sum(profile.w * d) for d in D)) / (2 * params.t0)
    phi = 2 * params.a * float(sum(np.sum(profile.w * (d - mean) ** 2) for d in D))
    return phi, 2 * mean

def cost(control, params, data, T, cells=256, breaks=()):
    profile = BaseProfile(params, data, T, cells=cells, breaks=breaks)
    total = 0.0
    for p in (0, 1):
        for j in profile.slots(p):
            total += float(np.sum(profile.w * _slot_values(control, profile, j) ** 2))
    return total
