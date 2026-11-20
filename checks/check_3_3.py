import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import sin, cos, pi
from profiles import make_profiles, fd_solve

def simpson(f, lo, hi, n=2000):
    h = (hi - lo) / n
    return h / 3 * (f(lo) + f(hi) + 4 * sum(f(lo + (2 * i - 1) * h) for i in range(1, n // 2 + 1))
                    + 2 * sum(f(lo + 2 * i * h) for i in range(1, n // 2)))

def family(l, a, y0, F):
    t0 = l / a
    def make(r, after=0.0):
        def u(t):
            if t < t0:
                return 0.5 * (y0(a * t) + F(a * t) / a + r)
            if t < 2 * t0:
                return 0.5 * (-y0(2 * l - a * t) + F(2 * l - a * t) / a + r)
            return after
        return u
    return make

def J_num(u, t0): # интеграл по двум половинам отдельно (в t0 возможен скачок)
    e = 1e-12
    return simpson(lambda t: u(t) ** 2, 0, t0 - e) + simpson(lambda t: u(t) ** 2, t0, 2 * t0 - e)


def spread(f, lo, hi, n=400):
    v = [f(lo + (hi - lo) * (i + 0.5) / n) for i in range(n)]
    return min(v), max(v)

l, a = 1.5, 2.0
t0 = l / a
T = 2 * t0
y0 = lambda x: 0.6 * (1 - x / l) + 0.7 * sin(pi * x / l) + 0.2 * sin(3 * pi * x / l)
y1 = lambda x: 0.5 * sin(2 * pi * x / l) + 0.3 * sin(pi * x / l)
F = lambda x: 0.5 * l / (2 * pi) * (1 - cos(2 * pi * x / l)) + 0.3 * l / pi * (1 - cos(pi * x / l))
make = family(l, a, y0, F)
xs = [l * i / 200 for i in range(1, 200)]

print("А) Формула сл. 7: y(2t0,x) = y0(x) + u(tau+) - u(tau-) для ПРОИЗВОЛЬНОГО u")
ub = lambda t: 0.6 * cos(3 * t) + 0.4 * sin(2 * t) ** 2
_, _, y = make_profiles(l, a, y0, F, ub)
print("   max расхождение = %.2e" % max(abs(y(T, x) - (y0(x) + ub(2 * t0 - x / a) - ub(x / a))) for x in xs))

print("Б) Семейство u_r: конечное состояние через профили")
for r in (y0(0), -0.37, 1.9, 0.0):
    A, B, y = make_profiles(l, a, y0, F, make(r))
    amin, amax = spread(A, T, T + t0)
    bmin, bmax = spread(B, T, T + t0)
    print("   r=%6.3f: max|y(T)|=%.1e, размах A=%.1e, размах B=%.1e, A=%.4f (= -r), B=%.4f"
          % (r, max(abs(y(T, x)) for x in xs), amax - amin, bmax - bmin, amin, bmin))

print("В) Независимая разностная схема (N=300 узлов)")
for r in (y0(0), -0.37, 1.9):
    dx, dt, n, prev, cur = fd_solve(l, a, y0, F, make(r), T)
    nxt = [cur[j + 1] + cur[j - 1] - prev[j] for j in range(1, len(cur) - 1)]
    vy = [abs(v) for v in cur]
    vt = [abs(nx - pv) / (2 * dt) for nx, pv in zip(nxt, prev[1:-1])]
    print("   r=%6.3f: max|y(T)| все узлы = %.1e, без 3 узлов у x=0: %.1e;  max|y_t(T)| без 3 узлов у x=0: %.1e"
          % (r, max(vy), max(vy[3:]), max(vt[3:])))

print("Г) Чем u_r НЕ является: добавка, не равная постоянной, нарушает покой")
r = 0.4
d1 = lambda t: make(r)(t) + 0.2
d2 = lambda t: make(r)(t) + 0.2 * sin(pi * t / t0)
for name, uu in (("u_r + 0.2          ", d1), ("u_r + 0.2 sin(pi t/t0)", d2)):
    _, _, y = make_profiles(l, a, y0, F, uu)
    A, B, _ = make_profiles(l, a, y0, F, uu)
    amin, amax = spread(A, T, T + t0)
    print("   %s: max|y(T)| = %.3e, размах A = %.3e" % (name, max(abs(y(T, x)) for x in xs), amax - amin))

print("Д) Стоимость: J(u_r) численно против формулы сл. 13 и разложения сл. 14")
Fbar = simpson(F, 0, l) / l
rstar = -Fbar / a
Jstar = simpson(lambda x: y0(x) ** 2 + (F(x) - Fbar) ** 2 / a ** 2, 0, l) / (2 * a)
print("   Fbar = %.6f, r* = %.6f, J* = %.8f" % (Fbar, rstar, Jstar))
for r in (rstar, y0(0), -0.37, 1.9):
    Jn = J_num(make(r), t0)
    J13 = simpson(lambda x: y0(x) ** 2 + (F(x) / a + r) ** 2, 0, l) / (2 * a)
    J14 = Jstar + l / (2 * a) * (r - rstar) ** 2
    print("   r=%9.6f: J числ.=%.8f, сл.13=%.8f, сл.14=%.8f" % (r, Jn, J13, J14))
eps = 1e-3
print("   минимум в r*: J(r*-e)-J(r*) = %.3e, J(r*+e)-J(r*) = %.3e (ожидается l/(2a) e^2 = %.3e)"
      % (J_num(make(rstar - eps), t0) - J_num(make(rstar), t0),
         J_num(make(rstar + eps), t0) - J_num(make(rstar), t0), l / (2 * a) * eps ** 2))

print("Е) Оптимальное h: интеграл по (0,2t0) и ортогональность")
h = make(rstar)
ih = simpson(h, 0, t0 - 1e-12) + simpson(h, t0, 2 * t0 - 1e-12)
print("   int h dt = %.2e" % ih)
r = 1.9
ur = make(r)
nd = simpson(lambda t: (ur(t) - h(t)) ** 2, 0, t0 - 1e-12) + simpson(lambda t: (ur(t) - h(t)) ** 2, t0, 2 * t0 - 1e-12)
print("   J(u_r)-J(h) = %.8f,  ||u_r-h||^2 = %.8f,  u_r-h = %.6f (ожидается (r-r*)/2 = %.6f)"
      % (J_num(ur, t0) - J_num(h, t0), nd, ur(0.3) - h(0.3), (r - rstar) / 2))

print("Ж) Согласование: r_en = y0(0)")
r = y0(0)
ur = make(r)
print("   u(0+) = %.6f, y0(0) = %.6f;  u(2t0-) = %.2e;  скачок в t0: %.2e (= -y0(l) = %.2e)"
      % (ur(1e-12), y0(0), ur(2 * t0 - 1e-12), ur(t0 + 1e-12) - ur(t0 - 1e-12), -y0(l)))

print("З) Пример лекции (сл. 19-20): y0 = kappa (l - x), y1 = 0")
kap = 0.8
z0 = lambda x: kap * (l - x)
zero = lambda x: 0.0
mk = family(l, a, z0, zero)
hh = mk(0.0)
print("   max|h - kappa/2 (l - a t)| = %.2e" % max(abs(hh(t) - kap / 2 * (l - a * t)) for t in [2 * t0 * i / 301 for i in range(1, 301)]))
Jh = J_num(hh, t0)
uen = mk(z0(0))
print("   J(h) = %.8f, kappa^2 l^3/(6a) = %.8f" % (Jh, kap ** 2 * l ** 3 / (6 * a)))
print("   max|u_en - (kappa l - kappa a t/2)| = %.2e,  J(u_en)/J(h) = %.6f"
      % (max(abs(uen(t) - (kap * l - kap * a * t / 2)) for t in [2 * t0 * i / 301 for i in range(1, 301)]), J_num(uen, t0) / Jh))
_, _, y = make_profiles(l, a, z0, zero, hh)
print("   max|y(2t0,x)| при h: %.2e" % max(abs(y(T, x)) for x in xs))

print("И) Свой пример: y0 = H sin(pi x/l), y1 = V")
H, V = 0.9, 1.3
w0 = lambda x: H * sin(pi * x / l)
FV = lambda x: V * x
mk = family(l, a, w0, FV)
rs = -V * t0 / 2
hh = mk(rs)
f1 = lambda t: 0.5 * (H * sin(pi * a * t / l) + V * (t - t0 / 2))
f2 = lambda t: 0.5 * (H * sin(pi * a * t / l) + V * (1.5 * t0 - t))
e = max(max(abs(hh(t) - f1(t)) for t in [t0 * i / 200 for i in range(1, 200)]),
        max(abs(hh(t) - f2(t)) for t in [t0 + t0 * i / 200 for i in range(1, 200)]))
print("   max|h - явная формула| = %.2e" % e)
print("   J(h) = %.8f, формула (1/2a)[H^2 l/2 + V^2 l^3/(12 a^2)] = %.8f"
      % (J_num(hh, t0), (H ** 2 * l / 2 + V ** 2 * l ** 3 / (12 * a ** 2)) / (2 * a)))
A, B, y = make_profiles(l, a, w0, FV, hh)
amin, amax = spread(A, T, T + t0)
print("   max|y(2t0,x)| = %.2e, размах A на (2t0,3t0) = %.2e" % (max(abs(y(T, x)) for x in xs), amax - amin))
