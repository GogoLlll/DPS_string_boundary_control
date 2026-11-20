import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import sin, cos, pi, floor, exp
from profiles import make_profiles

def simpson(f, lo, hi, n=400):
    if hi - lo < 1e-12 * max(1.0, abs(hi)):
        return 0.0
    d = (hi - lo) * 1e-11
    lo, hi = lo + d, hi - d
    h = (hi - lo) / n
    return h / 3 * (f(lo) + f(hi) + 4 * sum(f(lo + (2 * i - 1) * h) for i in range(1, n // 2 + 1))
                    + 2 * sum(f(lo + 2 * i * h) for i in range(1, n // 2)))

def integ(f, lo, hi, breaks):
    p = sorted({lo, hi} | {b for b in breaks if lo < b < hi})
    return sum(simpson(f, p[i], p[i + 1]) for i in range(len(p) - 1))

print("l a T/t0 max|y(T)|/масштаб J числ./J формула Phi(профили)/Phi(формула) пример лекции J/(k^2 l^3/6a)")
for l, a in ((1.0, 1.0), (2.3, 0.7), (0.4, 5.0), (7.0, 3.1)):
    t0 = l / a
    y0 = lambda x: 0.6 * (1 - x / l) + 0.7 * sin(pi * x / l) + 0.2 * sin(3 * pi * x / l) + 0.1 * exp(-x / l)
    y1 = lambda x: a / l * (0.5 * sin(2 * pi * x / l) + 0.3)
    F = lambda x: a / l * (0.5 * l / (2 * pi) * (1 - cos(2 * pi * x / l)) + 0.3 * x)
    def G(s):
        if s < t0:
            return 0.5 * (y0(a * s) + F(a * s) / a)
        return 0.5 * (-y0(2 * l - a * s) + F(2 * l - a * s) / a)
    n_of = lambda sig, T: 0 if sig >= T * (1 - 1e-13) else floor((T - sig) / (2 * t0) - 1e-12) + 1
    xs = [l * (i + 0.37) / 100 for i in range(0, 100)]
    kap = 0.8
    Gz = lambda s: kap / 2 * (l - a * s)
    Jz = simpson(lambda t: Gz(t) ** 2, 0, 2 * t0)
    for q in (0.8, 1.6, 2.0, 3.3, 6.0):
        T = q * t0
        rho = T - 2 * t0 * floor(T / (2 * t0) + 1e-12)
        br = [k * t0 for k in range(1, 20)] + [rho + 2 * k * t0 for k in range(0, 8)]
        if q >= 2:
            W = integ(lambda s: 1.0 / n_of(s, T), 0, 2 * t0, [t0, rho])
            M = integ(lambda s: G(s) / n_of(s, T), 0, 2 * t0, [t0, rho])
            S2 = integ(lambda s: G(s) ** 2 / n_of(s, T), 0, 2 * t0, [t0, rho])
            rs = -2 * M / W
            def u(t, T=T, rs=rs):
                if t < 0 or t >= T:
                    return 0.0
                sig = t - 2 * t0 * floor(t / (2 * t0))
                return (G(sig) + rs / 2) / n_of(sig, T)
            A, B, y = make_profiles(l, a, y0, F, u)
            Jn = integ(lambda t: u(t) ** 2, 0, T, br)
            print(" %5.1f  %5.1f   %.1f      %.1e            %.10f        -                          %.10f"
                  % (l, a, q, max(abs(y(T, x)) for x in xs), Jn / (S2 - M * M / W), Jz / (kap ** 2 * l ** 3 / (6 * a))))
        else:
            gm = integ(G, T, 2 * t0, [t0]) / (2 * t0 - T)
            pf = 2 * a * integ(lambda s: (G(s) - gm) ** 2, T, 2 * t0, [t0])
            u = lambda t, T=T, gm=gm: G(t) - gm if 0 <= t < T else 0.0
            A, B, y = make_profiles(l, a, y0, F, u)
            c = integ(lambda z: (A(z) - B(z)) / 2, T, T + t0, br) / t0
            pn = a / 2 * integ(lambda z: (A(z) - c) ** 2 + (B(z) + c) ** 2, T, T + t0, br)
            print(" %5.1f  %5.1f   %.1f      -                  -                   %.10f               -" % (l, a, q, pn / pf))
