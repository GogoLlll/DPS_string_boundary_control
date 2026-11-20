import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import sin, cos, pi, acos
import check_independent_sub as c

l, a, t0 = c.l, c.a, c.t0
c.CENTRAL = True

def basis(T, m):
    """Список функций: многочлены Чебышёва на каждом из кусков (0, min(T,t0)) и (t0, T)."""
    pieces = [(0.0, min(T, t0))] + ([(t0, T)] if T > t0 + 1e-12 else [])
    fs = []
    for lo, hi in pieces:
        for k in range(m + 1):
            def f(t, lo=lo, hi=hi, k=k):
                if not (lo <= t < hi or (hi == T and abs(t - T) < 1e-12)):
                    return 0.0
                return cos(k * acos(max(-1.0, min(1.0, 2 * (t - lo) / (hi - lo) - 1))))
            fs.append(f)
    return fs

def fit(q, m, N, y0, F):
    T = q * t0
    nT = round(q * N)
    dt = t0 / N
    fs = basis(T, m)
    base = c.residual_vector(N, nT, [0.0] * (nT + 1), y0, F)
    cols = [c.residual_vector(N, nT, [0.0] + [f(k * dt) for k in range(1, nT + 1)], None, None) for f in fs]
    rows, nc = len(base), len(fs)
    Cu = lambda v: [sum(cols[k][i] * v[k] for k in range(nc)) for i in range(rows)]
    Ct = lambda r: [sum(cols[k][i] * r[i] for i in range(rows)) for k in range(nc)]
    coef = c.cg(lambda v: Ct(Cu(v)), Ct([-v for v in base]), tol=1e-13, itmax=4000)
    res = [p + b for p, b in zip(Cu(coef), base)]
    u = lambda t: sum(ck * f(t) for ck, f in zip(coef, fs))
    return u, sum(v * v for v in res), sum(v * v for v in base)

if __name__ == "__main__":
    y0 = lambda x: 0.6 * (1 - x / l) + 0.7 * sin(pi * x / l) + 0.2 * sin(3 * pi * x / l)
    F = lambda x: 0.5 * l / (2 * pi) * (1 - cos(2 * pi * x / l)) + 0.3 * l / pi * (1 - cos(pi * x / l))
    kap = 0.8
    z0 = lambda x: kap * (l - x)
    zero = lambda x: 0.0
    mids = lambda lo, hi, n: [lo + (hi - lo) * (i + 0.5) / n for i in range(n)]
    for name, d0, dF in (("данные общего вида", y0, F), ("пример лекции y0 = kappa (l - x)", z0, zero)):
        print(name)
        print("  T/t0   степень  узлов   Phi без упр.   Phi_min найд.   Phi_min формула   max|u_найд - u_best|")
        for q in (0.5, 1.0, 1.5, 1.8):
            ub, pf = c.theory(q, d0, dF)
            for m, N in ((4, 200), (8, 200), (8, 400)):
                u, pmin, p0 = fit(q, m, N, d0, dF)
                T = q * t0
                pts = [t for t in mids(0, T, 300) if abs(t - t0) > 0.01 * t0]
                print("  %.1f     %d       %d     %.6f       %.6f        %.6f          %.2e"
                      % (q, m, N, p0, pmin, pf, max(abs(u(t) - ub(t)) for t in pts)))
