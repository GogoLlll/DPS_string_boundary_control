import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import sin, cos, pi, floor
from profiles import make_profiles, fd_solve

l, a = 1.5, 2.0
t0 = l / a
EPS = 1e-11
xs = [l * i / 200 for i in range(1, 200)]
mids = lambda lo, hi, n: [lo + (hi - lo) * (i + 0.5) / n for i in range(n)]

def simpson(f, lo, hi, n=400):
    if hi - lo < 1e-9:
        return 0.0
    lo, hi = lo + EPS, hi - EPS
    h = (hi - lo) / n
    return h / 3 * (f(lo) + f(hi) + 4 * sum(f(lo + (2 * i - 1) * h) for i in range(1, n // 2 + 1))
                    + 2 * sum(f(lo + 2 * i * h) for i in range(1, n // 2)))

def integ(f, lo, hi, breaks=()):
    p = sorted({lo, hi} | {b for b in breaks if lo < b < hi})
    return sum(simpson(f, p[i], p[i + 1]) for i in range(len(p) - 1))

def make_G(y0, F):
    def G(s):
        if s < t0:
            return 0.5 * (y0(a * s) + F(a * s) / a)
        return 0.5 * (-y0(2 * l - a * s) + F(2 * l - a * s) / a)
    return G

def phi_from_profiles(A, B, T, breaks):
    """Phi(T) = min_c (a/2) int_T^{T+t0} [(A-c)^2 + (B+c)^2] dz; минимум при c = среднее (A - B)/2."""
    br = [T + b for b in breaks] + [k * t0 for k in range(1, 8)]
    c = integ(lambda z: (A(z) - B(z)) / 2, T, T + t0, br) / t0
    return a / 2 * integ(lambda z: (A(z) - c) ** 2 + (B(z) + c) ** 2, T, T + t0, br), c

def phi_min_formula(G, T, breaks):
    if 2 * t0 - T < 1e-9:
        return 0.0, 0.0
    gm = integ(G, T, 2 * t0, breaks) / (2 * t0 - T)
    return 2 * a * integ(lambda s: (G(s) - gm) ** 2, T, 2 * t0, breaks), gm

print("А) Волна, бегущая влево, сосредоточена на (0, b): y0 = f, y1 = a f'")
b, H = 0.6 * l, 0.9
f = lambda x: H * cos(pi * x / (2 * b)) ** 2 if 0 <= x < b else 0.0
ya0 = f
Fa = lambda x: a * (f(x) - f(0))
Ga = make_G(ya0, Fa)
Tm = b / a
vals = [Ga(s) for s in mids(Tm, 2 * t0, 400)]
print("   T_min = b/a = %.3f t0;  G на (T_min, 2 t0): в [%.10f, %.10f], -y0(0)/2 = %.10f" % (Tm / t0, min(vals), max(vals), -ya0(0) / 2))
r = ya0(0)
u = lambda t: Ga(t) + r / 2 if 0 <= t < Tm else 0.0
print("   u = G + r/2 против f(a t): max расхождение = %.1e" % max(abs(u(t) - f(a * t)) for t in mids(0, Tm, 300)))
A, B, y = make_profiles(l, a, ya0, Fa, u)
print("   решение против f(x + a t): max расхождение = %.1e" % max(abs(y(t, x) - f(x + a * t)) for t in (0.1 * t0, 0.3 * t0, 0.5 * t0) for x in xs))
dx, dt, n, prev, cur = fd_solve(l, a, ya0, Fa, u, Tm)
nxt = [cur[j + 1] + cur[j - 1] - prev[j] for j in range(1, len(cur) - 1)]
print("   покой в T_min: профили max|y| = %.1e; разностная схема max|y| = %.1e, max|y_t| = %.1e"
      % (max(abs(y(Tm, x)) for x in xs), max(abs(v) for v in cur), max(abs(nx - pv) / (2 * dt) for nx, pv in zip(nxt, prev[1:-1]))))
for q in (0.9, 0.7):
    T = q * Tm
    print("   T = %.1f T_min: Phi_min = %.6f > 0 (точное успокоение невозможно)" % (q, phi_min_formula(Ga, T, [t0, Tm])[0]))

print("Б) Волна влево на всём (0, l) плюс волна вправо, сосредоточенная на (d, l); T_min = (2l - d)/a")
d = 0.45 * l
gL = lambda x: 0.5 * sin(pi * x / l) + 0.3 * sin(2 * pi * x / l) + 0.4 * (1 - x / l)       # левая, gL(l) = 0
dgL = lambda x: 0.5 * pi / l * cos(pi * x / l) + 0.6 * pi / l * cos(2 * pi * x / l) - 0.4 / l
psi = lambda x: 0.7 * sin(pi * (x - d) / (l - d)) ** 2 if d < x < l else 0.0               # правая
yb0 = lambda x: gL(x) + psi(x)
Fb = lambda x: a * (gL(x) - gL(0)) - a * psi(x)        # y1 = a gL' - a psi'
Gb = make_G(yb0, Fb)
Tm = (2 * l - d) / a
vals = [Gb(s) for s in mids(Tm, 2 * t0, 400)]
print("   T_min = %.3f t0;  G на (T_min, 2 t0): в [%.10f, %.10f], -y0(0)/2 = %.10f" % (Tm / t0, min(vals), max(vals), -yb0(0) / 2))
r = yb0(0)
u = lambda t: Gb(t) + r / 2 if 0 <= t < Tm else 0.0
A, B, y = make_profiles(l, a, yb0, Fb, u)
dx, dt, n, prev, cur = fd_solve(l, a, yb0, Fb, u, Tm, N=400)
nxt = [cur[j + 1] + cur[j - 1] - prev[j] for j in range(1, len(cur) - 1)]
print("   покой в T_min: профили max|y| = %.1e; разностная схема max|y| = %.1e, max|y_t| = %.1e"
      % (max(abs(y(Tm, x)) for x in xs), max(abs(v) for v in cur), max(abs(nx - pv) / (2 * dt) for nx, pv in zip(nxt, prev[1:-1]))))
print("   условие в координатах: max|y1 - a y0'| на (0, d) = %.1e;  на (d, l) = %.2f"
      % (max(abs((Fb(x + 1e-6) - Fb(x - 1e-6)) / 2e-6 - a * (yb0(x + 1e-6) - yb0(x - 1e-6)) / 2e-6) for x in mids(0.01, d - 0.01, 100)),
         max(abs((Fb(x + 1e-6) - Fb(x - 1e-6)) / 2e-6 - a * (yb0(x + 1e-6) - yb0(x - 1e-6)) / 2e-6) for x in mids(d, l, 100))))
for q in (1.4, 1.2, 0.8):
    print("   T = %.1f t0 < T_min: Phi_min = %.6f > 0" % (q, phi_min_formula(Gb, q * t0, [t0, Tm])[0]))

print("В) Тождество (a/2) int_T^{T+t0} [(A+r)^2 + (B-r)^2] dz = 2a int_0^{2t0} D_r^2 dsigma  (любые u, T, r)")
y0 = lambda x: 0.6 * (1 - x / l) + 0.7 * sin(pi * x / l) + 0.2 * sin(3 * pi * x / l)
F = lambda x: 0.5 * l / (2 * pi) * (1 - cos(2 * pi * x / l)) + 0.3 * l / pi * (1 - cos(pi * x / l))
G = make_G(y0, F)
w = lambda t: 0.6 * cos(3 * t) + 0.4 * sin(2 * t) ** 2 - 0.3 * sin(5 * t + 1)
def n_of(sig, T):
    return 0 if sig >= T else floor((T - sig) / (2 * t0) - 1e-12) + 1
for q, r in ((0.35, 0.3), (0.8, -0.5), (1.0, 0.0), (1.6, 0.7), (2.0, 0.2), (3.3, -0.4)):
    T = q * t0
    uw = lambda t: w(t) if 0 <= t < T else 0.0
    A, B, y = make_profiles(l, a, y0, F, uw)
    rho = T - 2 * t0 * floor(T / (2 * t0) + 1e-12)
    lhs = a / 2 * integ(lambda z: (A(z) + r) ** 2 + (B(z) - r) ** 2, T, T + t0, [k * t0 for k in range(1, 8)])
    D = lambda s: sum(w(s + 2 * m * t0) for m in range(n_of(s, T))) - G(s) - r / 2
    rhs = 2 * a * integ(lambda s: D(s) ** 2, 0, 2 * t0, [t0, rho])
    st = integ(lambda x: y(T, x) ** 2 + (0.5 * (A(T + x / a) - B(T + (l - x) / a)) + r) ** 2, 0, l,
               [l - a * (k * t0 - T) for k in range(1, 8)] + [a * (k * t0 - T) for k in range(1, 8)])
    print("   T = %.2f t0, r = %4.1f: через профили %.8f, через невязку %.8f, через состояние %.8f" % (q, r, lhs, rhs, st))

print("Г) Данные общего вида: управление u = G - Gср на (0, T) и его остаток")
print("   T/t0   Phi(u_best) числ.  формула      Phi(u=0)    Phi(u_best+возм.)  Phi(G + др. r)")
for q in (0.3, 0.7, 1.0, 1.3, 1.7, 1.95):
    T = q * t0
    pf, gm = phi_min_formula(G, T, [t0])
    ub = lambda t: G(t) - gm if 0 <= t < T else 0.0
    A, B, y = make_profiles(l, a, y0, F, ub)
    pn, c = phi_from_profiles(A, B, T, [])
    A0_, B0_, _ = make_profiles(l, a, y0, F, lambda t: 0.0)
    p0, _ = phi_from_profiles(A0_, B0_, T, [])
    up = lambda t: ub(t) + (0.2 * sin(pi * t / T) if 0 <= t < T else 0.0)
    Ap, Bp, _ = make_profiles(l, a, y0, F, up)
    pp, _ = phi_from_profiles(Ap, Bp, T, [])
    ur = lambda t: G(t) - gm + 0.25 if 0 <= t < T else 0.0
    Ar, Br, _ = make_profiles(l, a, y0, F, ur)
    pr, _ = phi_from_profiles(Ar, Br, T, [])
    print("   %.2f   %.8f         %.8f   %.6f    %.6f           %.6f" % (q, pn, pf, p0, pp, pr))

print("Д) Остаток как функция T")
vals = [phi_min_formula(G, 0.05 * i * t0, [t0])[0] for i in range(0, 41)]
Fbar = simpson(F, 0, l, 2000) / l
init = simpson(lambda x: y0(x) ** 2 + (F(x) - Fbar) ** 2 / a ** 2, 0, l, 2000)
print("   Phi_min(0) = %.8f, норма начального состояния int[y0^2 + (F-Fср)^2/a^2] = %.8f" % (vals[0], init))
print("   монотонно убывает: %s;  значения при T/t0 = 0, 0.5, 1, 1.5, 2: %s"
      % (all(vals[i + 1] <= vals[i] + 1e-12 for i in range(40)), ", ".join("%.5f" % vals[i] for i in (0, 10, 20, 30, 40))))

print("Е) Пример лекции y0 = kappa (l - x), y1 = 0: Phi_min(T) = kappa^2 a^3 (2 t0 - T)^3 / 24")
kap = 0.8
z0 = lambda x: kap * (l - x)
Gz = make_G(z0, lambda x: 0.0)
for q in (0.0, 0.5, 1.0, 1.5):
    T = q * t0
    pf, gm = phi_min_formula(Gz, T, [t0])
    line = "   T = %.1f t0: формула п. 3.6 %.8f, явная %.8f" % (q, pf, kap ** 2 * a ** 3 * (2 * t0 - T) ** 3 / 24)
    if q > 0:
        ub = lambda t: Gz(t) - gm if 0 <= t < T else 0.0
        A, B, y = make_profiles(l, a, z0, lambda x: 0.0, ub)
        line += ", через профили %.8f; u_best(t) - kappa a (t0 + T/2 - t)/2: %.1e" % (
            phi_from_profiles(A, B, T, [])[0], max(abs(ub(t) - kap * a * (t0 + T / 2 - t) / 2) for t in mids(0, T, 100)))
    print(line)
