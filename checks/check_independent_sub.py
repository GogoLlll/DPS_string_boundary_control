import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import sin, cos, pi, sqrt

l, a = 1.5, 2.0
t0 = l / a
CENTRAL = False

def residual_vector(N, nT, uvec, y0, F):
    dx = l / N
    dt = dx / a
    if y0 is not None:
        prev = [y0(j * dx) for j in range(N + 1)]
        cur = [0.0] * (N + 1)
        for j in range(1, N):
            cur[j] = 0.5 * (y0((j + 1) * dx) + y0((j - 1) * dx)) + (F((j + 1) * dx) - F((j - 1) * dx)) / (2 * a)
    else:
        prev = [0.0] * (N + 1)
        cur = [0.0] * (N + 1)
    cur[0] = uvec[1]
    for n in range(1, nT):
        nxt = [0.0] * (N + 1)
        for j in range(1, N):
            nxt[j] = cur[j + 1] + cur[j - 1] - prev[j]
        nxt[0] = uvec[n + 1]
        prev, cur = cur, nxt
    nxt = [0.0] * (N + 1)
    for j in range(1, N):
        nxt[j] = cur[j + 1] + cur[j - 1] - prev[j]
    if CENTRAL:
        yt = [(nxt[j] - prev[j]) / (2 * dt) for j in range(1, N)]
    else:
        yt = [(nxt[j] - cur[j]) / dt for j in range(1, N)]              # y_t(T, x_j), j = 1..N-1
    Q = [0.0]                                                          # Q в полуцелых точках x_{j+1/2}
    for v in yt:
        Q.append(Q[-1] + v * dx)
    qm = sum(Q) / len(Q)
    w = sqrt(dx)
    return [w * cur[j] for j in range(1, N)] + [w / a * (q - qm) for q in Q]

def cg(apply, b, tol=1e-14, itmax=20000):
    x = [0.0] * len(b)
    r = b[:]
    p = r[:]
    rr = sum(v * v for v in r)
    b2 = rr
    for _ in range(itmax):
        if rr <= tol * tol * max(b2, 1e-300):
            break
        Ap = apply(p)
        al = rr / sum(pi_ * q for pi_, q in zip(p, Ap))
        x = [xi + al * pi_ for xi, pi_ in zip(x, p)]
        r = [ri - al * q for ri, q in zip(r, Ap)]
        rr2 = sum(v * v for v in r)
        p = [ri + rr2 / rr * pi_ for ri, pi_ in zip(r, p)]
        rr = rr2
    return x

def minimize(N, q, y0, F):
    nT = round(q * N)
    base = residual_vector(N, nT, [0.0] * (nT + 1), y0, F)             # остаток без управления
    cols = []
    for k in range(1, nT + 1):
        e = [0.0] * (nT + 1)
        e[k] = 1.0
        cols.append(residual_vector(N, nT, e, None, None))             # линейная часть
    rows = len(base)
    Cu = lambda u: [sum(cols[k][i] * u[k] for k in range(nT)) for i in range(rows)]
    Ct = lambda r: [sum(cols[k][i] * r[i] for i in range(rows)) for k in range(nT)]
    u = cg(lambda v: Ct(Cu(v)), Ct([-v for v in base]))
    res = [p + b for p, b in zip(Cu(u), base)]
    return nT, u, sum(v * v for v in res), sum(v * v for v in base)

def simpson(f, lo, hi, n=400):
    if hi - lo < 1e-9:
        return 0.0
    lo, hi = lo + 1e-11, hi - 1e-11
    h = (hi - lo) / n
    return h / 3 * (f(lo) + f(hi) + 4 * sum(f(lo + (2 * i - 1) * h) for i in range(1, n // 2 + 1))
                    + 2 * sum(f(lo + 2 * i * h) for i in range(1, n // 2)))

def theory(q, y0, F):
    T = q * t0
    def G(s):
        if s < t0:
            return 0.5 * (y0(a * s) + F(a * s) / a)
        return 0.5 * (-y0(2 * l - a * s) + F(2 * l - a * s) / a)
    I = lambda f: simpson(f, T, t0) + simpson(f, max(T, t0), 2 * t0) if T < t0 else simpson(f, T, 2 * t0)
    gm = I(G) / (2 * t0 - T)
    return (lambda t: G(t) - gm), 2 * a * I(lambda s: (G(s) - gm) ** 2)

if __name__ == "__main__":
    y0 = lambda x: 0.6 * (1 - x / l) + 0.7 * sin(pi * x / l) + 0.2 * sin(3 * pi * x / l)
    F = lambda x: 0.5 * l / (2 * pi) * (1 - cos(2 * pi * x / l)) + 0.3 * l / pi * (1 - cos(pi * x / l))
    kap = 0.8
    z0 = lambda x: kap * (l - x)
    zero = lambda x: 0.0
    for name, d0, dF in (("данные общего вида", y0, F), ("пример лекции y0 = kappa (l - x)", z0, zero)):
        print(name)
        print("  T/t0    N   Phi без упр.  Phi_min дискр.  Phi_min формула   среднее |u_д - u_best|   пила (чёт - нечёт)")
        for q in (0.5, 1.0, 1.5, 1.8):
            ub, pf = theory(q, d0, dF)
            for N in (20, 40, 80):
                nT, u, pmin, p0 = minimize(N, q, d0, dF)
                dt = t0 / N
                d = [u[k - 1] - ub(k * dt) for k in range(1, nT - 3) if abs(k * dt - t0) > 1e-9]
                ev = [v for i, v in enumerate(d) if i % 2 == 0]
                od = [v for i, v in enumerate(d) if i % 2 == 1]
                print("  %.1f   %3d   %.6f      %.6f        %.6f         %.2e                %+.4f"
                      % (q, N, p0, pmin, pf, sum(abs(v) for v in d) / len(d), sum(ev) / len(ev) - sum(od) / len(od)))
    print("Для сравнения - центральная разность по времени (подрешётки не связаны), данные общего вида, T = 1.5 t0:")
    CENTRAL = True
    ub, pf = theory(1.5, y0, F)
    for N in (20, 40):
        nT, u, pmin, p0 = minimize(N, 1.5, y0, F)
        d = [u[k - 1] - ub(k * t0 / N) for k in range(1, nT - 3) if abs(k * t0 / N - t0) > 1e-9]
        ev = [v for i, v in enumerate(d) if i % 2 == 0]
        od = [v for i, v in enumerate(d) if i % 2 == 1]
        print("  N = %d: Phi_min = %.6f, среднее |u_д - u_best| = %.2e, пила = %+.4f"
              % (N, pmin, sum(abs(v) for v in d) / len(d), sum(ev) / len(ev) - sum(od) / len(od)))
