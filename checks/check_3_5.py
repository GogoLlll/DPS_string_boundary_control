import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import sin, cos, pi, floor, sqrt
from profiles import make_profiles, fd_solve

l, a = 1.5, 2.0
t0 = l / a
y0 = lambda x: 0.6 * (1 - x / l) + 0.7 * sin(pi * x / l) + 0.2 * sin(3 * pi * x / l)
y1 = lambda x: 0.5 * sin(2 * pi * x / l) + 0.3 * sin(pi * x / l)
F = lambda x: 0.5 * l / (2 * pi) * (1 - cos(2 * pi * x / l)) + 0.3 * l / pi * (1 - cos(pi * x / l))
EPS = 1e-11


def simpson(f, lo, hi, n=400):
    if hi - lo < 1e-9:
        return 0.0
    lo, hi = lo + EPS, hi - EPS
    h = (hi - lo) / n
    return h / 3 * (f(lo) + f(hi) + 4 * sum(f(lo + (2 * i - 1) * h) for i in range(1, n // 2 + 1))
                    + 2 * sum(f(lo + 2 * i * h) for i in range(1, n // 2)))


def integ(f, lo, hi, T):
    """Интеграл кусочно-гладкой функции: разрывы возможны в k t0 и в rho + 2 m t0."""
    rho = T - 2 * t0 * floor(T / (2 * t0) + 1e-12)
    pts = {lo, hi}
    k = 0
    while k * t0 < hi:
        if lo < k * t0 < hi: pts.add(k * t0)
        if lo < rho + k * t0 < hi: pts.add(rho + k * t0)   # с запасом
        k += 1
    p = sorted(pts)
    return sum(simpson(f, p[i], p[i + 1]) for i in range(len(p) - 1))


def make_G(y0, F):
    def G(s):
        if s < t0:
            return 0.5 * (y0(a * s) + F(a * s) / a)
        return 0.5 * (-y0(2 * l - a * s) + F(2 * l - a * s) / a)
    return G


def n_of(sig, T):
    return 0 if sig >= T else floor((T - sig) / (2 * t0) - 1e-12) + 1


def rstar_of(G, T):
    W = integ(lambda s: 1.0 / n_of(s, T), 0, 2 * t0, T)
    M = integ(lambda s: G(s) / n_of(s, T), 0, 2 * t0, T)
    S2 = integ(lambda s: G(s) ** 2 / n_of(s, T), 0, 2 * t0, T)
    return -2 * M / W, S2 - M * M / W, W, M


def min_norm(G, T, r):
    def u(t):
        if t < 0 or t >= T:
            return 0.0
        sig = t - 2 * t0 * floor(t / (2 * t0))
        return (G(sig) + r / 2) / n_of(sig, T)
    return u


G = make_G(y0, F)
xs = [l * i / 200 for i in range(1, 200)]
mids = lambda lo, hi, n: [lo + (hi - lo) * (i + 0.5) / n for i in range(n)]
J_of = lambda u, T: integ(lambda t: u(t) ** 2, 0, T, T)
Fbar = simpson(F, 0, l, 2000) / l
Jstar2 = simpson(lambda x: y0(x) ** 2 + (F(x) - Fbar) ** 2 / a ** 2, 0, l, 2000) / (2 * a)

print("А) Управление минимальной нормы: покой, стоимость, оптимальность по r")
print(" T/t0 r*(T) J числ. J формула  max|y(T)| размах A  J(r*-e)-J  J(r*+e)-J  (W/4)e^2")
for q in (2.0, 2.5, 3.0, 3.3, 4.0, 4.6, 5.7, 6.0, sqrt(31), 8.0):
    T = q * t0
    rs, Jf, W, M = rstar_of(G, T)
    u = min_norm(G, T, rs)
    A, B, y = make_profiles(l, a, y0, F, u)
    va = [A(z) for z in mids(T, T + t0, 200)]
    e = 1e-3
    Jn = J_of(u, T)
    print("   %5.3f  %9.6f  %.8f  %.8f   %.1e   %.1e   %.3e  %.3e  %.3e" % (
        q, rs, Jn, Jf, max(abs(y(T, x)) for x in xs), max(va) - min(va),
        J_of(min_norm(G, T, rs - e), T) - Jn, J_of(min_norm(G, T, rs + e), T) - Jn, W / 4 * e * e))

print("Б) T = 2 K t0: r* = -Fbar/a = %.6f, J = J*/K (J* = %.8f)" % (-Fbar / a, Jstar2))
for K in (1, 2, 3, 4):
    T = 2 * K * t0
    rs, Jf, W, M = rstar_of(G, T)
    u = min_norm(G, T, rs)
    h = min_norm(G, 2 * t0, -Fbar / a)
    e = max(abs(u(t) - h(t - 2 * t0 * floor(t / (2 * t0))) / K) for t in mids(0, T, 500))
    print("   K = %d: r* = %.6f, J = %.8f, J*/K = %.8f, max|u - h/K| = %.1e" % (K, rs, Jf, Jstar2 / K, e))

print("В) Минимальность при фиксированном r: J(u) = J(u_min) + ||u - u_min||^2 для других точных управлений")
T = 3.3 * t0
r = 0.6
um = min_norm(G, T, r)
v = lambda s: 0.5 * sin(pi * (s - 0.2 * t0) / (0.5 * t0)) ** 2 if 0.2 * t0 < s < 0.7 * t0 else 0.0
for name, uo in (("+v(t) - v(t-2t0)   ", lambda t: um(t) + v(t) - v(t - 2 * t0)),
                 ("+0.3 [1 на (0,rho), -1 на (2t0,2t0+rho)]", lambda t: um(t) + (0.3 if t < 1.3 * t0 else (-0.3 if 2 * t0 <= t < T else 0.0)))):
    A, B, y = make_profiles(l, a, y0, F, uo)
    d2 = integ(lambda t: (uo(t) - um(t)) ** 2, 0, T, T)
    print("   %s: max|y(T)| = %.1e, J(u) - J(u_min) = %.8f, ||u - u_min||^2 = %.8f"
          % (name, max(abs(y(T, x)) for x in xs), J_of(uo, T) - J_of(um, T), d2))

print("Г) Промежуточные состояния при T = 2 K t0: (y, y_t)(2 j t0) = (1 - j/K)(y0, y1)")
K = 4
T = 2 * K * t0
u = min_norm(G, T, rstar_of(G, T)[0])
A, B, y = make_profiles(l, a, y0, F, u)
hh = 1e-6
for j in range(K + 1):
    tj = 2 * j * t0
    ey = max(abs(y(tj, x) - (1 - j / K) * y0(x)) for x in xs)
    if j < K:
        ev = max(abs((y(tj + hh, x) - y(tj, x)) / hh - (1 - j / K) * y1(x)) for x in xs[5:-5])
    else:
        ev = max(abs((y(tj, x) - y(tj - hh, x)) / hh) for x in xs[5:-5])
    print("   j = %d: max|y - (1-j/K) y0| = %.1e, max|y_t - (1-j/K) y1| = %.1e (разностная производная)" % (j, ey, ev))

print("Д) Стоимость убывает с ростом T")
prev = None
ok = True
vals = []
for i in range(0, 61):
    q = 2.0 + 0.1 * i
    Jf = rstar_of(G, q * t0)[1]
    if prev is not None and Jf > prev + 1e-12: ok = False
    prev = Jf
    if i % 10 == 0: vals.append("%.1f: %.5f" % (q, Jf))
print("   J_min(T), T/t0 от 2 до 8 с шагом 0.1: монотонно не возрастает: %s;  %s" % (ok, "; ".join(vals)))

print("Е) Сравнение: 'погасить за 2 t0 оптимальным h, затем ноль' против минимальной нормы")
h = min_norm(G, 2 * t0, -Fbar / a)
for q in (3.0, 4.0, 6.0):
    T = q * t0
    uz = lambda t: h(t) if t < 2 * t0 else 0.0
    A, B, y = make_profiles(l, a, y0, F, uz)
    print("   T = %.1f t0: max|y(T)| = %.1e, J = %.6f против J_min = %.6f" % (q, max(abs(y(T, x)) for x in xs), Jstar2, rstar_of(G, T)[1]))

print("Ж) Непрерывный вариант: u_en (r = y0(0)) на (0, 2 t0), затем ноль")
uen = min_norm(G, 2 * t0, y0(0))
uz = lambda t: uen(t) if 0 <= t < 2 * t0 else 0.0
print("   u(0+) - y0(0) = %.1e, скачок в t0: %.1e, u(2t0-) = %.1e" % (uz(1e-9) - y0(0), uz(t0 + 1e-9) - uz(t0 - 1e-9), uz(2 * t0 - 1e-9)))
for q in (2.5, 3.3, 5.7):
    T = q * t0
    dx, dt, n, prev, cur = fd_solve(l, a, y0, F, uz, T)
    nxt = [cur[j + 1] + cur[j - 1] - prev[j] for j in range(1, len(cur) - 1)]
    print("   T = %.1f t0, разностная схема: max|y| = %.1e, max|y_t| = %.1e"
          % (q, max(abs(v) for v in cur), max(abs(nx - pv) / (2 * dt) for nx, pv in zip(nxt, prev[1:-1]))))
r = y0(0)
print("   стоимость J(u_en) = %.6f (не зависит от T)" % J_of(uen, 2 * t0))

print("З) Пример лекции y0 = kappa (l - x), y1 = 0 на горизонтах 3 t0 и 4 t0")
kap = 0.8
z0 = lambda x: kap * (l - x)
zero = lambda x: 0.0
Gz = make_G(z0, zero)
c = kap ** 2 * l ** 3 / a
T = 3 * t0
rs, Jf, W, M = rstar_of(Gz, T)
u = min_norm(Gz, T, rs)
f = lambda t: (kap / 4 * (l - a * (t if t < t0 else t - 2 * t0)) + kap * l / 24) if (t < t0 or t >= 2 * t0) else kap / 2 * (l - a * t) + kap * l / 12
A, B, y = make_profiles(l, a, z0, zero, u)
print("   T = 3 t0: r* = %.8f (kappa l/6 = %.8f), J = %.8f (11/96 kappa^2 l^3/a = %.8f), max|u - формула| = %.1e, max|y(T)| = %.1e"
      % (rs, kap * l / 6, Jf, 11 / 96 * c, max(abs(u(t) - f(t)) for t in mids(0, T, 600)), max(abs(y(T, x)) for x in xs)))
T = 4 * t0
rs, Jf, W, M = rstar_of(Gz, T)
u = min_norm(Gz, T, rs)
f = lambda t: kap / 4 * (l - a * (t - 2 * t0 * floor(t / (2 * t0))))
A, B, y = make_profiles(l, a, z0, zero, u)
print("   T = 4 t0: r* = %.1e, J = %.8f (kappa^2 l^3/(12a) = %.8f), max|u - формула| = %.1e, max|y(T)| = %.1e, max|y(2t0) - y0/2| = %.1e"
      % (rs, Jf, c / 12, max(abs(u(t) - f(t)) for t in mids(0, T, 600)), max(abs(y(T, x)) for x in xs),
         max(abs(y(2 * t0, x) - z0(x) / 2) for x in xs)))
