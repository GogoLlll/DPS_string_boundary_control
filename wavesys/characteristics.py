from dataclasses import dataclass

import numpy as np

from .segments import split_time

def first_layer_region(params, t, x):
    if not (0 <= t < params.t0):
        return None
    d_left, d_right = x / params.a, (params.l - x) / params.a
    if t < min(d_left, d_right):
        return "I"
    if d_left <= t < d_right:
        return "II"
    if d_right <= t < d_left:
        return "III"
    return "IV"

def point_info(params, t, x):
    t0 = params.t0
    zA, zB = t + x / params.a, t + (params.l - x) / params.a
    (kA, sA), (kB, sB) = split_time(zA, t0), split_time(zB, t0)
    kA, sA, kB, sB = int(kA), float(sA), int(kB), float(sB)
    times_A = [j * t0 + sA for j in range(kA % 2, kA, 2)]
    times_B = [j * t0 + sB for j in range((kB + 1) % 2, kB, 2)] 
    return {
        "zA": zA, "kA": kA, "sA": sA, "zB": zB, "kB": kB, "sB": sB,
        "control_times_A": times_A, "control_times_B": times_B,
        "region": first_layer_region(params, t, x),
        "A_from": "A0" if kA % 2 == 0 else "-B0", "B_from": "B0" if kB % 2 == 0 else "-A0",
    }

@dataclass
class Source:
    kind: str
    time: float
    size: float
    label: str = ""

def fold(distance, l):
    m = np.mod(distance, 2 * l)
    return np.where(m <= l, m, 2 * l - m)

def front_position(params, source, t):
    t = np.asarray(t, dtype=float)
    travelled = params.a * (t - source.time)
    start = 0.0 if source.kind == "left" else params.l
    pos = fold(start + travelled, params.l)
    return np.where(t >= source.time, pos, np.nan)

def front_polyline(params, source, horizon):
    t0 = params.t0
    first = source.time + (t0 if source.kind == "left" else t0)
    times = [source.time]
    k = 0
    while first + k * t0 < horizon:
        times.append(first + k * t0)
        k += 1
    times.append(horizon)
    times = np.array(times)
    return times, front_position(params, source, times)

def front_segments(params, sources, horizon, tol=1e-6):
    t0 = params.t0
    groups = {}
    for s in sources:
        time, size = (s.time, s.size) if s.kind == "left" else (s.time - t0, -s.size)
        phase = int(round((time % (2 * t0)) / (2 * t0) * 1e9)) % 10 ** 9
        groups.setdefault(phase, []).append((time, size))
    scale = max([1.0] + [abs(s.size) for s in sources])
    lines = []
    for events in groups.values():
        events.sort()
        total = 0.0
        for i, (time, size) in enumerate(events):
            total += size
            start = max(time, 0.0)
            end = events[i + 1][0] if i + 1 < len(events) else horizon
            end = min(end, horizon)
            if abs(total) <= tol * scale or end <= start:
                continue
            first = int(np.floor((start - time) / t0 + 1e-12)) + 1
            times = [start]
            k = first
            while time + k * t0 < end - 1e-12 * max(horizon, t0):
                times.append(time + k * t0)
                k += 1
            times.append(end)
            times = np.array(times)
            lines.append((times, front_position(params, Source("left", time, total), times), total))
    return lines

def find_sources(params, data, control, horizon, tol=1e-9):
    zero = np.array([0.0])
    u0 = float(control.slot(0, zero, params.t0)[0])
    y00 = float(np.asarray(data.y0(zero)).ravel()[0])
    y0l = float(np.asarray(data.y0(np.array([params.l]))).ravel()[0])
    scale = max(1.0, abs(u0), abs(y00), abs(y0l))
    sources = []
    if abs(u0 - y00) > tol * scale:
        sources.append(Source("left", 0.0, u0 - y00, "левый угол: u(0) ≠ y0(0)"))
    if abs(y0l) > tol * scale:
        sources.append(Source("right", 0.0, -y0l, "правый угол: y0(l) ≠ 0"))
    eps = 1e-9 * max(horizon, params.t0)
    for tau in control.jump_times(params.t0, horizon):
        size = float(control(np.array([tau + eps]), params.t0)[0] - control(np.array([tau - eps]), params.t0)[0])
        if abs(size) > tol * scale:
            sources.append(Source("left", float(tau), size, f"скачок управления при t = {tau:.6g}"))
    return sources
