// Вкладка «Движение»: форма струны в выбранный момент и её проигрывание во времени.
import { h, card } from "../dom.js";
import { fmt, inT0 } from "../format.js";
import { breakAt, frontPosition } from "../fronts.js";
import { colors, layout, line, plot } from "../theme.js";

let chart, slider, playBtn, speedSel, readout;
let data = null, frame = 0, playing = false, raf = 0, last = 0, acc = 0, lay = null, c = null;

function nearest(arr, v) {
  let best = 0;
  for (let i = 1; i < arr.length; i++) if (Math.abs(arr[i] - v) < Math.abs(arr[best] - v)) best = i;
  return best;
}

function draw() {
  if (!data) return;
  const p = data.params, t = data.grid.t[frame];
  const cuts = data.sources.map((s) => frontPosition(p, s, t));
  const shape = breakAt(data.grid.x, data.grid.Y[frame], cuts);
  const traces = [
    line(c, c.reference, { x: data.grid.x, y: data.y0, name: "начальное отклонение y₀(x)", line: { color: c.reference, width: 1.5 }, hoverinfo: "skip" }),
    line(c, c.s1, { x: shape.x, y: shape.y, name: "форма струны y(t, x)", hovertemplate: "y = %{y:.4g}<extra></extra>" }),
    { type: "scatter", mode: "markers", x: [0], y: [data.boundary.u[frame]], showlegend: false, hoverinfo: "skip",
      marker: { color: c.s1, size: 10, line: { color: c.surface, width: 2 } } },
  ];
  plot(chart, traces, lay);
  slider.value = String(frame);
  const iN = nearest(data.norm.t, t);
  readout.replaceChildren(
    h("span", {}, "t = ", h("b", {}, fmt(t)), ` (${inT0(t, p.t0)})`),
    h("span", {}, "управление u(t) = ", h("b", {}, fmt(data.boundary.u[frame]))),
    h("span", {}, "норма состояния Φ(t) = ", h("b", {}, fmt(data.norm.phi[iN]))),
  );
}

function tick(now) {
  if (!playing) return;
  const n = data.grid.t.length;
  const seconds = Math.min(60, Math.max(4, 2.5 * data.stats.T_over_t0)) / parseFloat(speedSel.value);
  acc += ((now - last) / 1000) * ((n - 1) / seconds);
  last = now;
  const step = Math.floor(acc);
  if (step > 0) {
    acc -= step;
    frame += step;
    if (frame >= n - 1) { frame = n - 1; draw(); stop(); return; }
    draw();
  }
  raf = requestAnimationFrame(tick);
}

function stop() {
  playing = false;
  cancelAnimationFrame(raf);
  playBtn.textContent = "▶ Пуск";
}

function toggle() {
  if (!data) return;
  if (playing) { stop(); return; }
  if (frame >= data.grid.t.length - 1) frame = 0;
  playing = true;
  playBtn.textContent = "❚❚ Пауза";
  last = performance.now();
  acc = 0;
  raf = requestAnimationFrame(tick);
}

export default {
  id: "motion",
  title: "Движение",
  mount(root) {
    chart = h("div", { class: "chart" });
    slider = h("input", { type: "range", min: "0", max: "0", value: "0", step: "1", "aria-label": "Момент времени",
      oninput: () => { stop(); frame = parseInt(slider.value, 10); draw(); } });
    playBtn = h("button", { type: "button", onclick: toggle }, "▶ Пуск");
    speedSel = h("select", { "aria-label": "Скорость" },
      ...[["0.25", "×0.25"], ["0.5", "×0.5"], ["1", "×1"], ["2", "×2"], ["4", "×4"]].map(([v, t]) => h("option", { value: v, selected: v === "1" }, t)));
    readout = h("div", { class: "readout" });
    root.append(card("Форма струны во времени",
      "Точка на левом конце — заданное управлением смещение. Там, где в решении скачок, линия прервана.",
      chart, h("div", { class: "controls" }, playBtn, h("label", {}, "скорость", speedSel), slider), readout));
  },
  update(d) {
    stop();
    data = d;
    c = colors();
    frame = 0;
    slider.max = String(d.grid.t.length - 1);
    const lo = Math.min(...d.grid.Y.map((row) => Math.min(...row)), ...d.y0);
    const hi = Math.max(...d.grid.Y.map((row) => Math.max(...row)), ...d.y0);
    const pad = 0.08 * Math.max(hi - lo, 1e-9);
    lay = layout(c, { x: "x", y: "y", legend: true, hovermode: "x unified",
      xaxis: { range: [0, d.params.l] }, yaxis: { range: [lo - pad, hi + pad] } });
    draw();
  },
  show() {
    // адрес с параметром play=1 сразу запускает проигрывание
    if (data && !playing && new URLSearchParams(window.location.hash.replace(/^#/, "")).get("play") === "1") toggle();
  },
  hide() { stop(); },
};
