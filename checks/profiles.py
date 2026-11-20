from math import floor

def make_profiles(l, a, y0, F, u):
    t0 = l / a
    A0 = lambda s: y0(a * s) + F(a * s) / a
    B0 = lambda s: y0(l - a * s) - F(l - a * s) / a

    def AB(k, s):
        Ak, Bk = A0(s), B0(s)
        for j in range(k):
            Ak, Bk = -Bk, 2 * u(j * t0 + s) - Ak
        return Ak, Bk

    def split(z):
        k = floor(z / t0)
        return k, z - k * t0

    A = lambda z: AB(*split(z))[0]
    B = lambda z: AB(*split(z))[1]
    y = lambda t, x: 0.5 * (A(t + x / a) + B(t + (l - x) / a))
    return A, B, y

def fd_solve(l, a, y0, F, u, T, N=300):
    dx = l / N
    dt = dx / a
    prev = [y0(j * dx) for j in range(N + 1)]
    cur = [0.0] * (N + 1)
    for j in range(1, N):
        cur[j] = 0.5 * (y0((j + 1) * dx) + y0((j - 1) * dx)) + (F((j + 1) * dx) - F((j - 1) * dx)) / (2 * a)
    cur[0] = u(dt)
    n = 1
    nT = round(T / dt)
    while n < nT:
        nxt = [0.0] * (N + 1)
        for j in range(1, N):
            nxt[j] = cur[j + 1] + cur[j - 1] - prev[j]
        n += 1
        nxt[0] = u(n * dt)
        prev, cur = cur, nxt
    return dx, dt, n, prev, cur
