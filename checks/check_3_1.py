from math import sin, cos, pi, floor
l, a = 1.5, 2.0
t0 = l / a
y0 = lambda x: 0.3*(1 - x/l) + 0.7*sin(pi*x/l) + 0.2*sin(3*pi*x/l)
y1 = lambda x: 0.5*sin(2*pi*x/l)
F  = lambda x: 0.5*l/(2*pi)*(1 - cos(2*pi*x/l))
u  = lambda t: 0.3*cos(3*t) + 0.4*sin(2*t)**2

A0 = lambda s: y0(a*s) + F(a*s)/a
B0 = lambda s: y0(l - a*s) - F(l - a*s)/a

def AB_rec(k, s):
    if k == 0: return A0(s), B0(s)
    Ap, Bp = AB_rec(k-1, s)
    return -Bp, 2*u((k-1)*t0 + s) - Ap

def AB_sum(k, s):
    m = k // 2
    if k % 2 == 0:
        return (A0(s) - 2*sum(u(s + 2*j*t0) for j in range(m)),
                B0(s) + 2*sum(u(s + (2*j+1)*t0) for j in range(m)))
    return (-B0(s) - 2*sum(u(s + (2*j+1)*t0) for j in range(m)),
            -A0(s) + 2*sum(u(s + 2*j*t0) for j in range(m+1)))

def split(z):
    k = floor(z / t0); return k, z - k*t0
A = lambda z: AB_rec(*split(z))[0]
B = lambda z: AB_rec(*split(z))[1]
y = lambda t, x: 0.5*(A(t + x/a) + B(t + (l - x)/a))

# 1) рекурсия == явные суммы
e = max(abs(p - q) for k in range(9) for i in range(1, 50)
        for p, q in zip(AB_rec(k, i*t0/50), AB_sum(k, i*t0/50)))
print("1) рекурсия vs явные суммы, k=0..8:        max|разн| = %.2e" % e)
# 2) начальные условия
h = 1e-6
xs = [l*i/200 for i in range(1, 200)]
print("2) y(0,x)-y0(x):                           max = %.2e" % max(abs(y(0, x) - y0(x)) for x in xs))
print("   y_t(0,x)-y1(x) (разность вперёд):       max = %.2e" % max(abs((y(h, x) - y(0, x))/h - y1(x)) for x in xs))
# 3) граничные условия до 6 t0
ts = [6*t0*i/997 for i in range(1, 997)]
print("3) y(t,0)-u(t), 0<t<6t0:                   max = %.2e" % max(abs(y(t, 0) - u(t)) for t in ts))
print("   y(t,l),      0<t<6t0:                   max = %.2e" % max(abs(y(t, l)) for t in ts))
# 4) уравнение: вторые разности во внутренних точках
d = 1e-3; e = 0
for t in [0.37*t0, 1.21*t0, 2.63*t0, 4.4*t0]:
    for x in [0.13*l, 0.5*l, 0.81*l]:
        ytt = (y(t+d, x) - 2*y(t, x) + y(t-d, x))/d**2
        yxx = (y(t, x+d) - 2*y(t, x) + y(t, x-d))/d**2
        e = max(e, abs(ytt - a*a*yxx))
print("4) y_tt - a^2 y_xx (конечные разности):     max = %.2e" % e)
# 5) независимое решение: схема "крест" с числом Куранта 1 (точна в узлах)
N = 300; dx = l/N; dt = dx/a
prev = [y0(j*dx) for j in range(N+1)]
cur = [0.0]*(N+1)
for j in range(1, N):   # первый слой по формуле Даламбера (точно)
    cur[j] = 0.5*(y0((j+1)*dx) + y0((j-1)*dx)) + (F((j+1)*dx) - F((j-1)*dx))/(2*a)
cur[0] = u(dt); cur[N] = 0.0
e = 0; n = 1
while n*dt < 5*t0:
    nxt = [0.0]*(N+1)
    for j in range(1, N): nxt[j] = cur[j+1] + cur[j-1] - prev[j]
    n += 1; nxt[0] = u(n*dt); nxt[N] = 0.0
    prev, cur = cur, nxt
    if n % 37 == 0:
        e = max(e, max(abs(cur[j] - y(n*dt, j*dx)) for j in range(1, N)))
print("5) профили vs разностная схема, 0<t<5t0:   max|разн| = %.2e" % e)
