import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import sin, cos, pi
from profiles import make_profiles, fd_solve

l, a = 1.5, 2.0
t0 = l / a
y0 = lambda x: 0.6 * (1 - x / l) + 0.7 * sin(pi * x / l) + 0.2 * sin(3 * pi * x / l)
y1 = lambda x: 0.5 * sin(2 * pi * x / l) + 0.3 * sin(pi * x / l)
F = lambda x: 0.5 * l / (2 * pi) * (1 - cos(2 * pi * x / l)) + 0.3 * l / pi * (1 - cos(pi * x / l))

def u_lecture6(r):
    def u(t):
        if t < t0:
            return 0.5 * (y0(a * t) + F(a * t) / a + r)
        if t < 2 * t0:
            return 0.5 * (-y0(2 * l - a * t) + F(2 * l - a * t) / a + r)
        return 0.0
    return u

def spread(f, lo, hi, n=400):
    v = [f(lo + (hi - lo) * (i + 0.5) / n) for i in range(n)]
    return min(v), max(v)

print("А) Точные управления лекции 6 (T = 2 t0): профили на (T, T+t0)")
for r in (y0(0), -0.37, 1.9):
    A, B, y = make_profiles(l, a, y0, F, u_lecture6(r))
    T = 2 * t0
    amin, amax = spread(A, T, T + t0)
    bmin, bmax = spread(B, T, T + t0)
    print("   r=%6.3f:  A в [%.12f, %.12f],  B в [%.12f, %.12f],  -r=%.3f" % (r, amin, amax, bmin, bmax, -r))

print("Б) То же управление, затем u = 0; момент T = 2.6 t0 (не кратен t0)")
r = y0(0)
u = u_lecture6(r)
A, B, y = make_profiles(l, a, y0, F, u)
T = 2.6 * t0
print("   A в [%.12f, %.12f],  B в [%.12f, %.12f]" % (spread(A, T, T + t0) + spread(B, T, T + t0)))

print("В) Независимая разностная схема при r = y0(0): состояние в T = 2 t0")
dx, dt, n, prev, cur = fd_solve(l, a, y0, F, u, 2 * t0)
nxt = [cur[j + 1] + cur[j - 1] - prev[j] for j in range(1, len(cur) - 1)]
print("   max|y(T,x)| = %.2e" % max(abs(v) for v in cur))
print("   max|y_t(T,x)| (центр. разность) = %.2e" % max(abs(nx - pv) / (2 * dt) for nx, pv in zip(nxt, prev[1:-1])))

print("Г) Обратная сторона: произвольное управление, профили не постоянны и состояние не нулевое")
ub = lambda t: 0.6 * cos(3 * t) + 0.4 * sin(2 * t) ** 2
A, B, y = make_profiles(l, a, y0, F, ub)
T = 2 * t0
amin, amax = spread(A, T, T + t0)
bmin, bmax = spread(B, T, T + t0)
ymax = max(abs(y(T, l * i / 200)) for i in range(1, 200))
print("   размах A = %.3f, размах B = %.3f, max|y(T,x)| = %.3f" % (amax - amin, bmax - bmin, ymax))

print("Д) Связь с характеристическими комбинациями: A' = y_t + a y_x, B' = y_t - a y_x")
h = 1e-5
e = 0
for t, x in [(0.4 * t0, 0.3 * l), (1.3 * t0, 0.55 * l), (2.2 * t0, 0.8 * l), (3.7 * t0, 0.15 * l)]:
    yt = (y(t + h, x) - y(t - h, x)) / (2 * h)
    yx = (y(t, x + h) - y(t, x - h)) / (2 * h)
    dA = (A(t + x / a + h) - A(t + x / a - h)) / (2 * h)
    dB = (B(t + (l - x) / a + h) - B(t + (l - x) / a - h)) / (2 * h)
    e = max(e, abs(dA - (yt + a * yx)), abs(dB - (yt - a * yx)))
print("   max расхождение = %.2e" % e)

print("Е) Сдвиг A -> A + c, B -> B - c (замена F на F + a c) не меняет решение")
c = 0.83
Fc = lambda x: F(x) + a * c
A2, B2, y2 = make_profiles(l, a, y0, Fc, ub)
pts = [(0.4 * t0, 0.3 * l), (1.3 * t0, 0.55 * l), (2.2 * t0, 0.8 * l), (3.7 * t0, 0.15 * l)]
print("   max|y - y2| = %.2e,  A2 - A = %.6f,  B2 - B = %.6f (ожидается +c и -c, c = %.2f)"
      % (max(abs(y(t, x) - y2(t, x)) for t, x in pts), A2(1.234) - A(1.234), B2(1.234) - B(1.234), c))
