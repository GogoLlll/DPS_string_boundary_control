import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import sin, cos, pi, floor
from profiles import make_profiles, fd_solve

l, a = 1.5, 2.0
t0 = l / a
EPS = 1e-11
xs = [l * i / 200 for i in range(1, 200)]
mids = lambda lo, hi, n: [lo + (hi - lo) * (i + 0.5) / n for i in range(n)]
y0 = lambda x: 0.6 * (1 - x / l) + 0.7 * sin(pi * x / l) + 0.2 * sin(3 * pi * x / l)
y1 = lambda x: 0.5 * sin(2 * pi * x / l) + 0.3 * sin(pi * x / l)
F = lambda x: 0.5 * l / (2 * pi) * (1 - cos(2 * pi * x / l)) + 0.3 * l / pi * (1 - cos(pi * x / l))

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

def G(s):
    if s < t0:
        return 0.5 * (y0(a * s) + F(a * s) / a)
    return 0.5 * (-y0(2 * l - a * s) + F(2 * l - a * s) / a)

def n_of(sig, T):
    return 0 if sig >= T else floor((T - sig) / (2 * t0) - 1e-12) + 1

def phi_link(u, T):
    """Способ 1: остаток по невязке условия связи (XXXIII), без решения прямой задачи."""
    rho = T - 2 * t0 * floor(T / (2 * t0) + 1e-12)
    D0 = lambda s: sum(u(s + 2 * m * t0) for m in range(n_of(s, T))) - G(s)
    mean = integ(D0, 0, 2 * t0, [t0, rho]) / (2 * t0)
    return 2 * a * integ(lambda s: (D0(s) - mean) ** 2, 0, 2 * t0, [t0, rho])

def phi_prof(A, B, t, breaks=()):
    """Способ 2: остаток через профили на (t, t + t0)."""
    br = list(breaks) + [k * t0 for k in range(0, 20)]
    c = integ(lambda z: (A(z) - B(z)) / 2, t, t + t0, br) / t0
    return a / 2 * integ(lambda z: (A(z) - c) ** 2 + (B(z) + c) ** 2, t, t + t0, br)

def phi_fd(u, T, N=600):
    """Способ 3: остаток по состоянию из разностной схемы: y и Q = int_0^x y_t dx (формула трапеций)."""
    dx, dt, n, prev, cur = fd_solve(l, a, y0, F, u, T, N=N)
    nxt = [0.0] * (N + 1)
    for j in range(1, N):
        nxt[j] = cur[j + 1] + cur[j - 1] - prev[j]
    nxt[0] = u(T + dt)
    yt = [(nxt[j] - prev[j]) / (2 * dt) for j in range(N + 1)]
    Q = [0.0]
    for j in range(1, N + 1):
        Q.append(Q[-1] + 0.5 * (yt[j] + yt[j - 1]) * dx)
    trap = lambda v: dx * (sum(v) - 0.5 * (v[0] + v[-1]))
    qm = trap(Q) / l
    return trap([v * v for v in cur]) + trap([(q - qm) ** 2 for q in Q]) / a ** 2

print("А) Три независимых способа вычислить остаток Phi(T) для гладкого НЕуспокаивающего управления")
w = lambda t: 0.6 * cos(3 * t) + 0.4 * sin(2 * t) ** 2          # w(0) = y0(0), w'(0) = y1(0)
A, B, y = make_profiles(l, a, y0, F, w)
print("   T/t0    по невязке связи   через профили     разностная схема (N=600)")
for q in (0.5, 1.3, 2.0, 2.7, 4.5):
    T = q * t0
    print("   %.1f     %.8f         %.8f        %.8f" % (q, phi_link(w, T), phi_prof(A, B, T), phi_fd(w, T)))

print("Б) Тождество ромба y(t+dt,x) + y(t-dt,x) = y(t,x+dx) + y(t,x-dx), dx = a dt, для решения через профили")
import random
random.seed(1)
for name, uu in (("гладкое управление", w),
                 ("разрывное управление", lambda t: (G(t - 2 * t0 * floor(t / (2 * t0))) - 0.1) / 2 if 0 <= t < 4 * t0 else 0.0)):
    Ar, Br, yr = make_profiles(l, a, y0, F, uu)
    e = 0.0
    for _ in range(4000):
        dt_ = random.uniform(0.001, 0.2) * t0
        dx_ = a * dt_
        x = random.uniform(dx_, l - dx_)
        t = random.uniform(dt_, 5 * t0)
        e = max(e, abs(yr(t + dt_, x) + yr(t - dt_, x) - yr(t, x + dx_) - yr(t, x - dx_)))
    print("   %s, 4000 случайных ромбов: max невязка = %.1e" % (name, e))

print("В) Кривая Phi(t) вдоль траектории: управление минимальной нормы, T = 2 K t0, K = 3")
K = 3
T = 2 * K * t0
Fbar = simpson(F, 0, l, 2000) / l
u = lambda t: (G(t - 2 * t0 * floor(t / (2 * t0))) - Fbar / (2 * a)) / K if 0 <= t < T else 0.0
A, B, y = make_profiles(l, a, y0, F, u)
init = simpson(lambda x: y0(x) ** 2 + (F(x) - Fbar) ** 2 / a ** 2, 0, l, 2000)
print("   Phi(0) = %.8f, норма начального состояния = %.8f" % (phi_prof(A, B, 0.0), init))
for j in range(K + 1):
    print("   t = %d t0: Phi = %.8f, (1 - j/K)^2 Phi(0) = %.8f" % (2 * j, phi_prof(A, B, 2 * j * t0), (1 - j / K) ** 2 * init))
print("   в промежуточные моменты: " + ", ".join("t=%.1f t0: %.4f" % (q, phi_prof(A, B, q * t0)) for q in (0.5, 1.0, 1.5, 2.5, 3.5, 5.0, 5.5)))

print("Г) Покой сохраняется после T при u = 0")
print("   max|y(t,x)| при T < t < T + 3 t0: %.1e;  Phi(T + 1.37 t0) = %.1e"
      % (max(abs(y(T + q * t0, x)) for q in (0.01, 0.4, 1.0, 1.7, 2.9) for x in xs), phi_prof(A, B, T + 1.37 * t0)))

print("Д) Причинность: изменение управления после момента tau не меняет y(t, x) при t < tau")
tau = 1.7 * t0
u2 = lambda t: u(t) + (0.8 * sin(7 * t) if t >= tau else 0.0)
A2, B2, y2 = make_profiles(l, a, y0, F, u2)
print("   max|y - y2| при t < tau: %.1e;  при t = tau + 0.5 t0: %.3f"
      % (max(abs(y(t, x) - y2(t, x)) for t in mids(0, tau, 40) for x in xs),
         max(abs(y(tau + 0.5 * t0, x) - y2(tau + 0.5 * t0, x)) for x in xs)))

print("Е) Сводная проверка всех режимов тремя способами (остаток должен быть нулевым)")
cases = []
Tc = 2 * t0
cases.append(("T = 2 t0, h          ", Tc, lambda t: G(t) - Fbar / (2 * a) if 0 <= t < Tc else 0.0))
cases.append(("T = 2 t0, u_en       ", Tc, lambda t: G(t) + y0(0) / 2 if 0 <= t < Tc else 0.0))
cases.append(("T = 6 t0, u* = h/3   ", T, u))
Tz = 3.5 * t0
cases.append(("T = 3.5 t0, u_en, 0  ", Tz, lambda t: G(t) + y0(0) / 2 if 0 <= t < Tc else 0.0))
for name, TT, uu in cases:
    Ax, Bx, yx = make_profiles(l, a, y0, F, uu)
    dx, dt, n, prev, cur = fd_solve(l, a, y0, F, uu, TT)
    nbad = sum(1 for v in cur if abs(v) > 1e-9)
    print("   %s: по связи %.1e, через профили %.1e, max|y(T)| профили %.1e, схема: узлов с |y|>1e-9: %d из %d"
          % (name, phi_link(uu, TT), phi_prof(Ax, Bx, TT), max(abs(yx(TT, x)) for x in xs), nbad, len(cur)))
