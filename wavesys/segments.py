import numpy as np

SNAP = 1e-12

def split_time(z, t0):
    z = np.asarray(z, dtype=float)
    k = np.floor(z / t0).astype(np.int64)
    s = z - k * t0
    neg = s < 0
    k = np.where(neg, k - 1, k)
    s = np.where(neg, s + t0, s)
    over = s >= t0
    k = np.where(over, k + 1, k)
    s = np.where(over, s - t0, s)
    near_hi = s > t0 * (1 - SNAP)
    k = np.where(near_hi, k + 1, k)
    s = np.where(near_hi, 0.0, s)
    s = np.where(s < t0 * SNAP, 0.0, s)
    k = np.maximum(k, 0)
    return k, s


def before(j, s, k_end, s_end):
    return (j < k_end) | ((j == k_end) & (s < s_end))
