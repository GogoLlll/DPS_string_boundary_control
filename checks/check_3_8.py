import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import sin, cos, pi, floor
from profiles import make_profiles, fd_solve

l, a = 1.5, 2.0
t0 = l / a
E = 1e-9
xs = [l * i / 200 for i in range(1, 200)]
mids = lambda lo, hi, n: [lo + (hi - lo) * (i + 0.5) / n for i in range(n)]
y0 = lambda x: 0.6 * (1 - x / l) + 0.7 * sin(pi * x / l) + 0.2 * sin(3 * pi * x / l)
dy0 = lambda x: -0.6 / l + 0.7 * pi / l * cos(pi * x / l) + 0.6 * pi / l * cos(3 * pi * x / l)
y1 = lambda x: 0.5 * sin(2 * pi * x / l) + 0.3 * sin(pi * x / l)
F = lambda x: 0.5 * l / (2 * pi) * (1 - cos(2 * pi * x / l)) + 0.3 * l / pi * (1 - cos(pi * x / l))


def simpson(f, lo, hi, n=400):
    if hi - lo < 1e-9:
        return 0.0
    lo, hi = lo + 1e-11, hi - 1e-11
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


G = make_G(y0, F)
jump = lambda f, z: f(z + E) - f(z - E)

print("А) Закон для скачков профилей в точках z = k t0:  a_{k+1} = -b_k,  b_{k+1} = 2[u]_{k t0} - a_k")
# управление со скачками в моменты, кратные t0, и несогласованное в углу; y0(l) != 0
y0n = lambda x: y0(x) + 0.35 * x / l                      # y0n(l) = 0.35
vals = [0.9, -0.2, 0.5, 0.1, -0.6, 0.3]
un = lambda t: (vals[int(floor(t / t0))] + 0.3 * sin(2 * t)) if 0 <= t < 6 * t0 else 0.0
A, B, y = make_profiles(l, a, y0n, F, un)
ak, bk = -2 * y0n(l), 2 * (un(E) - y0n(0))
e = max(abs(jump(A, t0) - ak), abs(jump(B, t0) - bk))
for k in range(1, 6):
    ak, bk = -bk, 2 * jump(un, k * t0) - ak
    e = max(e, abs(jump(A, (k + 1) * t0) - ak), abs(jump(B, (k + 1) * t0) - bk))
print("   a_1 = -2 y0(l), b_1 = 2 (u(0+) - y0(0)); сравнение предсказанных и измеренных скачков для k = 1..6: max расхождение = %.1e" % e)

print("Б) Скачок решения вдоль характеристики: семейство u_r при T = 2 t0, r != y0(0)")
for r in (-0.37, 1.9):
    u = lambda t, r=r: G(t) + r / 2 if 0 <= t < 2 * t0 else 0.0
    A, B, y = make_profiles(l, a, y0, F, u)
    d0 = 0.5 * (r - y0(0))
    x = 0.37 * l
    j1 = y(x / a + E, x) - y(x / a - E, x)                       # прямой фронт t = x/a
    j2 = y((2 * l - x) / a + E, x) - y((2 * l - x) / a - E, x)   # отражённый фронт t = (2l - x)/a
    # вне фронтов решение непрерывно по времени
    cont = max(abs(y(t + E, x) - y(t - E, x)) for t in mids(0.01 * t0, 1.99 * t0, 400)
               if abs(t - x / a) > 1e-3 and abs(t - (2 * l - x) / a) > 1e-3)
    print("   r = %5.2f: delta0 = (r - y0(0))/2 = %.4f; скачок на t = x/a: %.4f; на t = (2l-x)/a: %.4f; вне фронтов: %.1e; u(2t0-) = %.4f"
          % (r, d0, j1, j2, cont, u(2 * t0 - E)))

print("В) u_en (r = y0(0)): решение непрерывно, но на характеристике t = x/a есть излом")
u = lambda t: G(t) + y0(0) / 2 if 0 <= t < 2 * t0 else 0.0
A, B, y = make_profiles(l, a, y0, F, u)
x = 0.37 * l
cont = max(abs(y(t + E, x) - y(t - E, x)) for t in mids(0.001 * t0, 2.5 * t0, 2000))
h = 1e-5
tf = x / a
yt_after = (y(tf + 2 * h, x) - y(tf + h, x)) / h
yt_before = (y(tf - h, x) - y(tf - 2 * h, x)) / h
print("   max скачок y по времени в точке x = 0.37 l: %.1e" % cont)
print("   скачок y_t на t = x/a: измерено %.5f; формула (a y0'(0) - y1(0))/2 = %.5f" % (yt_after - yt_before, 0.5 * (a * dy0(0) - y1(0))))

print("Г) Правый угол: y0(l) != 0. Скачок -y0(l) идёт от (0,l) влево и гасится скачком u_r в момент t0")
Gn = make_G(y0n, F)
u = lambda t: Gn(t) + y0n(0) / 2 if 0 <= t < 2 * t0 else 0.0
A, B, y = make_profiles(l, a, y0n, F, u)
x = 0.37 * l
tl = (l - x) / a
print("   скачок u в t0: %.4f (= -y0(l) = %.4f); скачок y на t = (l-x)/a: %.4f; скачок B в 2 t0: %.1e; max|y(2t0,x)| = %.1e"
      % (jump(u, t0), -y0n(l), y(tl + E, x) - y(tl - E, x), jump(B, 2 * t0), max(abs(y(2 * t0, xx)) for xx in xs)))

print("Д) Непрерывные управления при T = 2 K t0: веса omega_0 = 1 - (K-1) phi/K, omega_m = phi/K, phi(sigma) = min(sigma/eps, 1)")
g = lambda s: G(s) + y0(0) / 2                              # r = r_en: g(0+) = y0(0), g(2t0-) = 0
Juen = integ(lambda t: g(t) ** 2, 0, 2 * t0, [t0])
Fbar = simpson(F, 0, l, 2000) / l
Jstar = simpson(lambda x: y0(x) ** 2 + (F(x) - Fbar) ** 2 / a ** 2, 0, l, 2000) / (2 * a)
for K in (2, 3):
    T = 2 * K * t0
    print("   K = %d: J_min = J*/K = %.5f;  J(u_en)/K = %.5f;  J(u_en, затем ноль) = %.5f" % (K, Jstar / K, Juen / K, Juen))
    for eps in (1.0, 0.3, 0.05):
        ee = eps * t0
        phi = lambda s: min(s / ee, 1.0)
        def u(t, K=K, T=T, phi=phi):
            if t < 0 or t >= T:
                return 0.0
            m = int(floor(t / (2 * t0)))
            s = t - 2 * m * t0
            w = 1 - (K - 1) * phi(s) / K if m == 0 else phi(s) / K
            return w * g(s)
        maxjump = max(abs(u(t + 1e-7) - u(t - 1e-7)) for t in [T * i / 4000 for i in range(1, 4000)] + [k * t0 for k in range(1, 2 * K)])
        J = integ(lambda t: u(t) ** 2, 0, T, [k * t0 for k in range(1, 2 * K)] + [2 * m * t0 + ee for m in range(K)])
        A, B, y = make_profiles(l, a, y0, F, u)
        dx, dt, n, prev, cur = fd_solve(l, a, y0, F, u, T)
        nxt = [cur[j + 1] + cur[j - 1] - prev[j] for j in range(1, len(cur) - 1)]
        print("      eps = %.2f t0: u(0+) - y0(0) = %.1e, max скачок u = %.1e, J = %.5f, покой: профили %.1e, схема |y| %.1e, |y_t| %.1e"
              % (eps, u(1e-12) - y0(0), maxjump, J, max(abs(y(T, x)) for x in xs), max(abs(v) for v in cur),
                 max(abs(nx - pv) / (2 * dt) for nx, pv in zip(nxt, prev[1:-1]))))

print("Е) Скачки управления минимальной нормы u* при T = 3.3 t0 (K = 1, rho = 1.3 t0)")
def n_of(sig, T):
    return 0 if sig >= T else floor((T - sig) / (2 * t0) - 1e-12) + 1
T = 3.3 * t0
rho = 1.3 * t0
W = integ(lambda s: 1.0 / n_of(s, T), 0, 2 * t0, [t0, rho])
M = integ(lambda s: G(s) / n_of(s, T), 0, 2 * t0, [t0, rho])
rs = -2 * M / W
us = lambda t: (G(t - 2 * t0 * floor(t / (2 * t0))) + rs / 2) / n_of(t - 2 * t0 * floor(t / (2 * t0)), T) if 0 <= t < T else 0.0
print("   r* = %.5f; u*(0+) - y0(0) = %.4f; скачки: t0: %.1e, rho: %.4f (формула (G(rho)+r*/2)/2 = %.4f), 2 t0: %.4f, 3 t0: %.1e, u*(T-) = %.4f"
      % (rs, us(E) - y0(0), jump(us, t0), jump(us, rho), (G(rho) + rs / 2) / 2, jump(us, 2 * t0), jump(us, 3 * t0), us(T - E)))
