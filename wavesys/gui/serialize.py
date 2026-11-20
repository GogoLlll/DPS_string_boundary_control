import math

import numpy as np

def num(v):
    if v is None:
        return None
    v = float(v)
    return v if math.isfinite(v) else None


def arr(a):
    a = np.asarray(a, dtype=float)
    if a.ndim == 0:
        return num(a)
    out = a.astype(object)
    out[~np.isfinite(a)] = None
    return out.tolist()
