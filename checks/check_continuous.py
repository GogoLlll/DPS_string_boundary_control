import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import sin, cos, pi, floor
from profiles import make_profiles, fd_solve

l, a = 1.5, 2.0
t0 = l / a
y0 = lambda x: 0.6 * (1 - x / l) + 0.7 * sin(pi * x / l) + 0.2 * sin(3 * pi * x / l)     # y0(l) = 0
F = lambda x: 0.5 * l / (2 * pi) * (1 - cos(2 * pi * x / l)) + 0.3 * l / pi * (1 - cos(pi * x / l))
xs = [l * (i + 0.37) / 100 for i in range(100)]

def G(s):
    if s < t0:
        return 0.5 * (y0(a * s) + F(a * s) / a)
    return 0.5 * (-y0(2 * l - a * s) + F(2 * l - a * s) / a)

def simpson(f, lo, hi, n=200):
    if hi - lo < 1e-10:
        return 0.0
    lo, hi = lo + 1e-12, hi - 1e-12
    h = (hi - lo) / n
    return h / 3 * (f(lo) + f(hi) + 4 * sum(f(lo + (2 * i - 1) * h) for i in range(1, n // 2 + 1))
                    + 2 * sum(f(lo + 2 * i * h) for i in range(1, n // 2)))

def integ(f, lo, hi, breaks):
    p = sorted({lo, hi} | {b for b in breaks if lo + 1e-10 < b < hi - 1e-10})
    return sum(simpson(f, p[i], p[i + 1]) for i in range(len(p) - 1))

def split_T(T):
    K = int(floor(T / (2 * t0) + 1e-9))
    rho = T - 2 * K * t0
    if rho < 1e-9:
        rho = 0.0
    return K, rho

def rstar_Jmin(T):
    K, rho = split_T(T)
    n = lambda s: K + 1 if s < rho else K
    W = integ(lambda s: 1.0 / n(s), 0, 2 * t0, [t0, rho])
    M = integ(lambda s: G(s) / n(s), 0, 2 * t0, [t0, rho])
    S2 = integ(lambda s: G(s) ** 2 / n(s), 0, 2 * t0, [t0, rho])
    return -2 * M / W, S2 - M * M / W, (lambda r: S2 + r * M + r * r * W / 4)

def make_cont(T, r, eps):
    """Непрерывное управление (XLII). Возвращает u(t) и список точек излома для интегрирования."""
    K, rho = split_T(T)
    g = lambda s: G(s) + r / 2
    g0 = 0.5 * (y0(0) + r)                  # g(0+)
    g2 = 0.5 * (-y0(0) + r)                 # g(2 t0 -)
    n0 = K + 1 if rho > 0 else K            # число повторов у sigma = 0+
    Nj = K if rho > 0 else K - 1            # число внутренних стыков t = 2 m t0
    assert Nj >= 1, "нужно T > 2 t0"
    vj = (g0 - y0(0)) / Nj
    def v(m):                               # значение управления в момент 2 m t0
        if m == 0:
            return y0(0)
        if rho == 0 and m == K:
            return 0.0
        return vj
    e = min(eps, rho, 2 * t0 - rho) if rho > 0 else min(eps, t0)
    grho = g(rho) if rho > 0 else 0.0
    def u(t):
        if t < 0 or t >= T:
            return 0.0
        m = int(floor(t / (2 * t0) + 1e-12))
        s = t - 2 * m * t0
        if s < 0:
            s = 0.0
        val = g(s) / (K + 1 if s < rho else K)
        if s < e:                                           # слой у sigma = 0+
            val += (v(m) - g0 / n0) * (1 - s / e)
        if 2 * t0 - s < e and m <= K - 1:                   # слой у sigma = 2 t0 -
            val += (v(m + 1) - g2 / K) * (1 - (2 * t0 - s) / e)
        if rho > 0 and s < rho and rho - s < e:             # слой у sigma = rho -
            wm = grho / K if m < K else 0.0
            val += (wm - grho / (K + 1)) * (1 - (rho - s) / e)
        return val
    br = [k * t0 for k in range(0, 2 * K + 3)]
    for m in range(K + 1):
        br += [2 * m * t0 + e, 2 * (m + 1) * t0 - e, 2 * m * t0 + rho, 2 * m * t0 + rho - e]
    return u, br

def n_of(sig, T):
    return 0 if sig >= T else int(floor((T - sig) / (2 * t0) - 1e-12)) + 1

if __name__ == "__main__":
    print("Непрерывные управления при r = r*(T): стоимость стремится к J_min(T) при eps -> 0")
    print("  T/t0   eps/t0   u(0+)-y0(0)   u(T-)     макс. скачок u   невязка связи   профили |y(T)|   схема |y|, |y_t|     J          J_min")
    for q in (2.5, 3.3, 4.0, 5.7, 6.0):
        T = q * t0
        rs, Jmin, Jr = rstar_Jmin(T)
        for eps in (0.3, 0.1, 0.02):
            u, br = make_cont(T, rs, eps * t0)
            step = T / 20000
            mj = max(abs(u((i + 1) * step) - u(i * step)) for i in range(1, 19998))
            mj = max(mj, max(abs(u(b + 1e-9) - u(b - 1e-9)) for b in br if 1e-6 < b < T - 1e-6))
            link = max(abs(sum(u(s + 2 * m * t0) for m in range(n_of(s, T))) - G(s) - rs / 2)
                       for s in [2 * t0 * (i + 0.5) / 400 for i in range(400)])
            A, B, y = make_profiles(l, a, y0, F, u)
            dx, dt, n, prev, cur = fd_solve(l, a, y0, F, u, T)
            nxt = [cur[j + 1] + cur[j - 1] - prev[j] for j in range(1, len(cur) - 1)]
            J = integ(lambda t: u(t) ** 2, 0, T, br)
            print("  %.1f    %.2f     %+.1e      %+.1e   %.1e          %.1e         %.1e          %.1e, %.1e     %.5f    %.5f"
                  % (q, eps, u(1e-10) - y0(0), u(T - 1e-10), mj, link, max(abs(y(T, x)) for x in xs),
                     max(abs(v_) for v_ in cur), max(abs(nx - pv) / (2 * dt) for nx, pv in zip(nxt, prev[1:-1])), J, Jmin))
    print("Для сравнения при T = 4 t0: J_min = %.5f; J(u_en)/2 = %.5f (к нему стремилось прежнее семейство с r = y0(0))"
          % (rstar_Jmin(4 * t0)[1], rstar_Jmin(4 * t0)[2](y0(0))))
    print("Произвольное r (не оптимальное), T = 3.3 t0: непрерывность и покой")
    T = 3.3 * t0
    for r in (-0.8, 0.6, 2.0):
        u, br = make_cont(T, r, 0.1 * t0)
        A, B, y = make_profiles(l, a, y0, F, u)
        dx, dt, n, prev, cur = fd_solve(l, a, y0, F, u, T)
        print("  r = %+.1f: u(0+)-y0(0) = %+.1e, u(T-) = %+.1e, профили |y(T)| = %.1e, схема |y| = %.1e, J = %.5f, формула int g^2/n = %.5f"
              % (r, u(1e-10) - y0(0), u(T - 1e-10), max(abs(y(T, x)) for x in xs), max(abs(v_) for v_ in cur),
                 integ(lambda t: u(t) ** 2, 0, T, br), rstar_Jmin(T)[2](r)))
