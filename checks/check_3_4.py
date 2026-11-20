import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import sin, cos, pi, floor, sqrt
from profiles import make_profiles, fd_solve

l, a = 1.5, 2.0
t0 = l / a
y0 = lambda x: 0.6 * (1 - x / l) + 0.7 * sin(pi * x / l) + 0.2 * sin(3 * pi * x / l)
y1 = lambda x: 0.5 * sin(2 * pi * x / l) + 0.3 * sin(pi * x / l)
F = lambda x: 0.5 * l / (2 * pi) * (1 - cos(2 * pi * x / l)) + 0.3 * l / pi * (1 - cos(pi * x / l))

def G(s):
    if s < t0:
        return 0.5 * (y0(a * s) + F(a * s) / a)
    return 0.5 * (-y0(2 * l - a * s) + F(2 * l - a * s) / a)

def n_of(sig, T):
    return 0 if sig >= T else floor((T - sig) / (2 * t0) - 1e-12) + 1

def stacked(u, sig, T):
    return sum(u(sig + 2 * m * t0) for m in range(n_of(sig, T)))

w = lambda t: 0.6 * cos(3 * t) + 0.4 * sin(2 * t) ** 2 - 0.3 * sin(5 * t + 1)   # произвольная функция
mids = lambda lo, hi, n: [lo + (hi - lo) * (i + 0.5) / n for i in range(n)]

print("А) Тождество для ПРОИЗВОЛЬНОГО u и любого T: на (T, T+t0)")
print("   A = -2 D(.), B = +2 D(.), где D(sigma) = sum u(sigma+2m t0) - G(sigma)")
A, B, y = make_profiles(l, a, y0, F, w)
for q in (0.35, 0.7, 1.0, 1.4, 2.0, 2.5, 3.0, 3.3, 4.0, 5.7, sqrt(31)):
    T = q * t0
    D = lambda sig: stacked(w, sig, T) - G(sig)
    e = 0.0
    for z in mids(T, T + t0, 300):
        k = floor(z / t0)
        s = z - k * t0
        sa, sb = (s, s + t0) if k % 2 == 0 else (s + t0, s)
        e = max(e, abs(A(z) + 2 * D(sa)), abs(B(z) - 2 * D(sb)))
    print("   T = %6.3f t0: max расхождение = %.2e" % (q, e))

def exact_control(T, r, filler):
    def u(t):
        if t >= T or t < 0:
            return 0.0
        m = floor(t / (2 * t0))
        sig = t - 2 * m * t0
        n = n_of(sig, T)
        if m < n - 1:
            return filler(t)
        return G(sig) + r / 2 - sum(filler(sig + 2 * j * t0) for j in range(n - 1))
    return u

print("Б) Достаточность: управление, построенное по условию связи (остальные слоты произвольны)")
xs = [l * i / 200 for i in range(1, 200)]
for q, r in ((2.0, 0.6), (2.5, 0.6), (3.0, -0.4), (3.3, 0.6), (4.0, 1.1), (5.7, 0.6), (sqrt(31), -0.2)):
    T = q * t0
    u = exact_control(T, r, w)
    A, B, y = make_profiles(l, a, y0, F, u)
    va = [A(z) for z in mids(T, T + t0, 300)]
    vb = [B(z) for z in mids(T, T + t0, 300)]
    line = "   T = %6.3f t0, r = %4.1f: max|y(T)| = %.1e, A в [%.6f, %.6f], B в [%.6f, %.6f]" % (
        q, r, max(abs(y(T, x)) for x in xs), min(va), max(va), min(vb), max(vb))
    if abs(q * 10 - round(q * 10)) < 1e-9:
        dx, dt, n, prev, cur = fd_solve(l, a, y0, F, u, T)
        big = sorted(abs(v) for v in cur)
        nbad = sum(1 for v in big if v > 1e-9)
        line += " | схема: узлов с |y| > 1e-9: %d из %d, в остальных max|y| = %.1e" % (
            nbad, len(cur), big[len(big) - 1 - nbad])
    print(line)

print("В) Необходимость: нарушение условия связи в одном слоте нарушает покой")
T = 3.3 * t0
base = exact_control(T, 0.6, w)
bump = lambda t: base(t) + (0.3 * sin(pi * (t - 0.2 * t0) / (0.5 * t0)) ** 2 if 0.2 * t0 < t < 0.7 * t0 else 0.0)
A, B, y = make_profiles(l, a, y0, F, bump)
print("   добавка на (0.2 t0, 0.7 t0): max|y(T)| = %.3f" % max(abs(y(T, x)) for x in xs))
v = lambda sig: 0.5 * sin(pi * (sig - 0.2 * t0) / (0.5 * t0)) ** 2 if 0.2 * t0 < sig < 0.7 * t0 else 0.0
ker = lambda t: base(t) + v(t) - v(t - 2 * t0)
A, B, y = make_profiles(l, a, y0, F, ker)
print("   добавка +v(t) - v(t - 2 t0) (суммы не меняются): max|y(T)| = %.1e,  max|y(1.5 t0)| отличие траектории = %.3f"
      % (max(abs(y(T, x)) for x in xs),
         max(abs(y(1.5 * t0, x) - make_profiles(l, a, y0, F, base)[2](1.5 * t0, x)) for x in xs)))

print("Г) Второй вывод через один профиль P (п. 2.9)")
Gt = lambda x: F(l) - F(x)
def P0(s):
    if s < 0:
        return 0.5 * (-y0(l + a * s) + Gt(l + a * s) / a)
    return 0.5 * (y0(l - a * s) + Gt(l - a * s) / a)
e = max(abs(-P0(sig - t0) - (G(sig) - F(l) / (2 * a))) for sig in mids(0, 2 * t0, 400))
print("   -P0(sigma - t0) = G(sigma) - F(l)/(2a): max расхождение = %.2e" % e)
def P(z, u):                                       # (2.33)
    m = floor((z + t0) / (2 * t0))
    s = z - 2 * m * t0
    return P0(s) + sum(u(s + (2 * j + 1) * t0) for j in range(m))
A, B, y = make_profiles(l, a, y0, F, w)
e = 0.0
for t, x in [(0.4 * t0, 0.3 * l), (1.3 * t0, 0.55 * l), (2.2 * t0, 0.8 * l), (3.7 * t0, 0.15 * l), (5.1 * t0, 0.62 * l)]:
    e = max(e, abs(y(t, x) - (P(t + (l - x) / a, w) - P(t - (l - x) / a, w))))
print("представление (2.29) против (2.5): max расхождение = %.2e" % e)
T = 3.3 * t0
u = exact_control(T, 0.6, w)
vp = [P(z, u) for z in mids(T - t0, T + t0, 400)]
print("   P на (T - t0, T + t0) при точном управлении, T = 3.3 t0: значения в [%.8f, %.8f]" % (min(vp), max(vp)))

print("Д) Частные случаи вручную")
T = 3 * t0
u = exact_control(T, 0.0, w)
s = 0.37 * t0
A0 = y0(a * s) + F(a * s) / a
B0 = y0(l - a * s) - F(l - a * s) / a
print("   T = 3 t0: u0+u2 - A0/2 = %.2e,  u1 + B0/2 = %.2e" % (u(s) + u(2 * t0 + s) - A0 / 2, u(t0 + s) + B0 / 2))

print("Е) Непрерывное управление (нет скачков) - чистая проверка разностной схемой при T, не кратном 2 t0")
y0c = lambda x: 0.7 * sin(pi * x / l) + 0.2 * sin(3 * pi * x / l)
def Gc(s):
    s = s % (2 * t0)
    if s < t0:
        return 0.5 * (y0c(a * s) + F(a * s) / a)
    return 0.5 * (-y0c(2 * l - a * s) + F(2 * l - a * s) / a)
def make_smooth(T):
    L = T - 2 * t0
    def Hc(tau):                                   # гладкая ступенька от 0 до 1 на [0, L]
        if tau <= 0: return 0.0
        if tau >= L: return 1.0
        return tau / L - sin(2 * pi * tau / L) / (2 * pi)
    om = lambda t: Hc(t) - Hc(t - 2 * t0)          # веса: сумма по сдвигам на 2 t0 равна 1
    return (lambda t: om(t) * Gc(t) if 0 <= t < T else 0.0), om
for q in (2.5, 3.3, 4.0, 5.7):
    T = q * t0
    u, om = make_smooth(T)
    es = max(abs(sum(u(sig + 2 * m * t0) for m in range(n_of(sig, T))) - Gc(sig)) for sig in mids(0, 2 * t0, 300))
    A, B, y = make_profiles(l, a, y0c, F, u)
    dx, dt, n, prev, cur = fd_solve(l, a, y0c, F, u, T)
    nxt = [cur[j + 1] + cur[j - 1] - prev[j] for j in range(1, len(cur) - 1)]
    print("   T = %.1f t0: невязка условия связи %.1e; профили: max|y(T)| = %.1e; схема: max|y| = %.1e, max|y_t| = %.1e"
          % (q, es, max(abs(y(T, x)) for x in xs), max(abs(v) for v in cur),
             max(abs(nx - pv) / (2 * dt) for nx, pv in zip(nxt, prev[1:-1]))))
