import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import sin, cos, pi, floor

l, a = 1.5, 2.0
t0 = l / a
y0 = lambda x: 0.6 * (1 - x / l) + 0.7 * sin(pi * x / l) + 0.2 * sin(3 * pi * x / l)
F = lambda x: 0.5 * l / (2 * pi) * (1 - cos(2 * pi * x / l)) + 0.3 * l / pi * (1 - cos(pi * x / l))

def run(N, nT, uvec, with_data, strict=False):
    """uvec[n] = u(t_n), n = 1..nT. strict=False: слои nT и nT+1; strict=True: u_nT = 0, слои nT-1 и nT."""
    dx = l / N
    if with_data:
        prev = [y0(j * dx) for j in range(N + 1)]
        cur = [0.0] * (N + 1)
        for j in range(1, N):
            cur[j] = 0.5 * (y0((j + 1) * dx) + y0((j - 1) * dx)) + (F((j + 1) * dx) - F((j - 1) * dx)) / (2 * a)
    else:
        prev = [0.0] * (N + 1)
        cur = [0.0] * (N + 1)
    cur[0] = uvec[1]
    last = nT - 1 if strict else nT
    for n in range(1, last + 1):
        nxt = [0.0] * (N + 1)
        for j in range(1, N):
            nxt[j] = cur[j + 1] + cur[j - 1] - prev[j]
        nxt[0] = uvec[n + 1] if n + 1 <= (nT - 1 if strict else nT) else 0.0
        prev, cur = cur, nxt
    return prev[1:N] + cur[1:N]

def cg(apply, b, tol=1e-13, itmax=5000):
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

def solve(N, q, strict=False):
    nT = round(q * N)                       # T = q t0, dt = t0 / N
    m = nT - 1 if strict else nT            # число неизвестных
    cols = []
    for k in range(1, m + 1):
        e = [0.0] * (nT + 2)
        e[k] = 1.0
        cols.append(run(N, nT, e, False, strict))   # отклик на единичное управление в момент t_k
    b = [-v for v in run(N, nT, [0.0] * (nT + 2), True, strict)]
    rows = len(b)
    Mu = lambda u: [sum(cols[k][i] * u[k] for k in range(m)) for i in range(rows)]
    Mt = lambda lam: [sum(cols[k][i] * lam[i] for i in range(rows)) for k in range(m)]
    # минимальная норма при M u = b:  u = M^T lam,  (M M^T) lam = b
    lam = cg(lambda v: Mu(Mt(v)), b)
    u = Mt(lam)
    res = max(abs(p - q_) for p, q_ in zip(Mu(u), b))
    return nT, u, res, rows, m


def G(s):
    if s < t0:
        return 0.5 * (y0(a * s) + F(a * s) / a)
    return 0.5 * (-y0(2 * l - a * s) + F(2 * l - a * s) / a)


def n_of(sig, T):
    return 0 if sig >= T else floor((T - sig) / (2 * t0) - 1e-12) + 1


def simpson(f, lo, hi, n=400):
    if hi - lo < 1e-9:
        return 0.0
    lo, hi = lo + 1e-11, hi - 1e-11
    h = (hi - lo) / n
    return h / 3 * (f(lo) + f(hi) + 4 * sum(f(lo + (2 * i - 1) * h) for i in range(1, n // 2 + 1))
                    + 2 * sum(f(lo + 2 * i * h) for i in range(1, n // 2)))


def theory(T):
    rho = T - 2 * t0 * floor(T / (2 * t0) + 1e-12)
    br = sorted({0.0, t0, 2 * t0} | ({rho} if 0 < rho < 2 * t0 else set()))
    I = lambda f: sum(simpson(f, br[i], br[i + 1]) for i in range(len(br) - 1))
    W = I(lambda s: 1.0 / n_of(s, T))
    M = I(lambda s: G(s) / n_of(s, T))
    S2 = I(lambda s: G(s) ** 2 / n_of(s, T))
    rs = -2 * M / W
    def u(t):
        sig = t - 2 * t0 * floor(t / (2 * t0) + 1e-12)
        return (G(sig) + rs / 2) / n_of(sig, T)
    return u, S2 - M * M / W, rs


if __name__ == "__main__":
    print("Дискретная задача на минимум нормы (без профилей и без условия связи) против формул (XV), (XXIV)")
    print("  T/t0    N   уравн. неизв.  невязка   J дискр.   J_min теория  max|u_д - u*|  среднее    u_д(T) и u*(T-)")
    for q in (2.0, 2.5, 3.0, 4.0, 5.0):
        for N in (20, 40, 80):
            if abs(q * N - round(q * N)) > 1e-9 or q * N > 250:
                continue
            nT, u, res, rows, m = solve(N, q)
            dt = t0 / N
            T = q * t0
            uth, Jth, rs = theory(T)
            J = sum(v * v for v in u) * dt
            rho = T - 2 * t0 * floor(T / (2 * t0) + 1e-12)
            brk = [k * t0 for k in range(0, 12)] + [rho + 2 * k * t0 for k in range(0, 6)]
            diffs = [abs(u[k - 1] - uth(k * dt)) for k in range(1, nT) if all(abs(k * dt - c) > 1e-9 for c in brk)]
            print("  %.1f   %3d   %4d  %4d   %.1e   %.6f   %.6f     %.2e      %.2e   %.4f  %.4f"
                  % (q, N, rows, m, res, J, Jth, max(diffs), sum(diffs) / len(diffs), u[-1], uth(T - 1e-9)))
    print("Для сравнения - постановка с запретом u(T) = 0 ('strict') при T = 2 t0:")
    for N in (20, 40):
        nT, u, res, rows, m = solve(N, 2.0, strict=True)
        print("  N = %d: J = %.6f (теория без запрета: %.6f)" % (N, sum(v * v for v in u) * t0 / N, theory(2 * t0)[1]))

    # ---------- докритический случай: разрешима ли дискретная система ----------
    print("Докритический случай: наименьшая невязка дискретной системы покоя (метод наименьших квадратов)")
    from math import cos as _cos
    b_, H_ = 0.6 * l, 0.9
    f_ = lambda x: H_ * _cos(pi * x / (2 * b_)) ** 2 if 0 <= x < b_ else 0.0
    g = globals()
    saved = (g["y0"], g["F"])

    def lsq(N, q):
        nT = round(q * N)
        cols = []
        for k in range(1, nT + 1):
            e = [0.0] * (nT + 2)
            e[k] = 1.0
            cols.append(run(N, nT, e, False))
        b = [-v for v in run(N, nT, [0.0] * (nT + 2), True)]
        rows, m = len(b), nT
        Mu = lambda u: [sum(cols[k][i] * u[k] for k in range(m)) for i in range(rows)]
        Mt = lambda lam: [sum(cols[k][i] * lam[i] for i in range(rows)) for k in range(m)]
        u = cg(lambda v: Mt(Mu(v)), Mt(b))            # нормальные уравнения M^T M u = M^T b
        r = [p - q_ for p, q_ in zip(Mu(u), b)]
        return u, (sum(v * v for v in r) * l / N) ** 0.5

    g["y0"], g["F"] = f_, (lambda x: a * (f_(x) - f_(0)))       # волна влево на (0, 0.6 l): T_min = 0.6 t0
    for q in (0.6, 0.5, 0.4):
        for N in (20, 40):
            u, rn = lsq(N, q)
            line = "  волна влево, T = %.1f t0, N = %d: невязка = %.2e" % (q, N, rn)
            if q == 0.6:
                line += ";  max|u_д - f(a t)| = %.2e" % max(abs(u[k - 1] - f_(a * k * t0 / N)) for k in range(1, round(q * N) + 1))
            print(line)
    g["y0"], g["F"] = saved
    for q in (1.0, 1.5, 2.0):
        u, rn = lsq(40, q)
        print("  данные общего вида, T = %.1f t0, N = 40: невязка = %.2e" % (q, rn))
