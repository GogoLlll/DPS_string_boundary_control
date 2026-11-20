// Линии разрыва решения. Нужны только для рисования: кривую нельзя соединять через скачок.
//
// Источник: { kind: "left" | "right", time, size }. Фронт выходит из левого (правого) конца в момент time,
// бежит со скоростью a и отражается от концов.

export function frontPosition(p, src, t) {
  if (t < src.time) return null;
  const start = src.kind === "left" ? 0 : p.l;
  const m = (start + p.a * (t - src.time)) % (2 * p.l);
  return m <= p.l ? m : 2 * p.l - m;
}

// Моменты из (0, T], когда фронт проходит через точку x.
export function crossingTimes(p, src, x, T) {
  const out = [];
  const first = src.kind === "left" ? [x, 2 * p.l - x] : [p.l - x, p.l + x];
  for (const d of first) {
    for (let m = 0; ; m++) {
      const t = src.time + (d + 2 * p.l * m) / p.a;
      if (t > T) break;
      out.push(t);
    }
  }
  return out;
}

// Аргументы z, в которых скачок имеют профили A и B (закон распространения скачков).
export function profileJumps(p, src, zMax) {
  const A = [], B = [];
  const firstB = src.kind === "left" ? src.time + p.t0 : src.time + 2 * p.t0;
  const firstA = src.kind === "left" ? src.time + 2 * p.t0 : src.time + p.t0;
  const end = zMax * (1 + 1e-12);
  for (let z = firstB; z <= end; z += 2 * p.t0) B.push(z);
  for (let z = firstA; z <= end; z += 2 * p.t0) A.push(z);
  return { A, B };
}

// Разрывает ломаную в заданных позициях.
//  - Если позиция лежит между соседними точками, обе части доводятся до неё значениями соседних точек
//    (приближение односторонних пределов), и между ними вставляется пустое значение.
//  - Если позиция совпадает с точкой сетки, эта точка выбрасывается: значение ровно на линии разрыва
//    относится к одной из сторон по соглашению и на графике только мешает.
export function breakAt(xs, ys, positions) {
  const lo = xs[0], hi = xs[xs.length - 1];
  const eps = 1e-9 * (hi - lo);
  const cuts = positions.filter((v) => v !== null && v >= lo - eps && v <= hi + eps).sort((a, b) => a - b);
  if (!cuts.length) return { x: xs, y: ys };
  const x = [], y = [];
  let c = 0, lastY = null, pending = null;
  for (let i = 0; i < xs.length; i++) {
    while (c < cuts.length && cuts[c] < xs[i] - eps) {          // разрыв между точками i-1 и i
      if (i > 0) {
        if (lastY !== null) { x.push(cuts[c]); y.push(lastY); }
        x.push(cuts[c]); y.push(null);
        pending = cuts[c];
      }
      c++;
    }
    if (c < cuts.length && Math.abs(cuts[c] - xs[i]) <= eps) {  // точка лежит на разрыве
      if (lastY !== null) { x.push(xs[i]); y.push(lastY); }
      x.push(xs[i]); y.push(null);
      while (c < cuts.length && Math.abs(cuts[c] - xs[i]) <= eps) c++;
      lastY = null; pending = xs[i];
      continue;
    }
    if (pending !== null && ys[i] !== null) { x.push(pending); y.push(ys[i]); }
    pending = null;
    x.push(xs[i]); y.push(ys[i]);
    lastY = ys[i];
  }
  return { x, y };
}
