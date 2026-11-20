// Вкладка «Сечения»: форма струны в несколько моментов и движение нескольких точек.
import { getJSON } from "../api.js";
import { h, card } from "../dom.js";
import { fmt, inT0 } from "../format.js";
import { breakAt, crossingTimes, frontPosition } from "../fronts.js";
import { colors, layout, line, ordinalColors, plot } from "../theme.js";

let shapesChart, tracesChart, timesInput, pointsInput, message, data = null;

async function load(ctx) {
  if (!data) return;
  message.textContent = "";
  try {
    const res = await getJSON("/api/sections", { times: timesInput.value, points: pointsInput.value });
    const c = colors(), p = data.params;
    const col1 = ordinalColors(c, res.shapes.length);
    plot(shapesChart, res.shapes.map((s, i) => {
      const cut = breakAt(res.x, s.y, data.sources.map((src) => frontPosition(p, src, s.t)));
      return line(c, col1[i], { x: cut.x, y: cut.y, name: `t = ${fmt(s.t)} (${inT0(s.t, p.t0)})`, hovertemplate: "%{y:.4g}" });
    }), layout(c, { x: "x", y: "y", legend: res.shapes.length > 1, xaxis: { range: [0, p.l] } }));
    const col2 = ordinalColors(c, res.traces.length);
    plot(tracesChart, res.traces.map((s, i) => {
      const cut = breakAt(res.t, s.y, data.sources.flatMap((src) => crossingTimes(p, src, s.x, p.T)));
      return line(c, col2[i], { x: cut.x, y: cut.y, name: `x = ${fmt(s.x)} (${fmt(s.x / p.l)} l)`, hovertemplate: "%{y:.4g}" });
    }), layout(c, { x: "t", y: "y", legend: res.traces.length > 1, xaxis: { range: [0, p.T] } }));
  } catch (err) {
    message.textContent = err.message;
  }
}

export default {
  id: "sections",
  title: "Сечения",
  mount(root, ctx) {
    timesInput = h("input", { type: "text", size: "34", "aria-label": "Моменты времени" });
    pointsInput = h("input", { type: "text", size: "34", "aria-label": "Точки струны" });
    message = h("p", { class: "hint", role: "alert" });
    shapesChart = h("div", { class: "chart" });
    tracesChart = h("div", { class: "chart" });
    const go = h("button", { type: "button", onclick: () => load(ctx) }, "Построить");
    root.append(
      h("div", { class: "card" },
        h("div", { class: "controls" },
          h("label", {}, "моменты t:", timesInput), h("label", {}, "точки x:", pointsInput), go),
        h("p", { class: "hint" }, "До пяти значений через запятую. Более тёмная линия — более поздний момент или более далёкая от левого конца точка."),
        message),
      card("Форма струны y(t, x) в заданные моменты", null, shapesChart),
      card("Движение заданных точек: y(t, x) как функция времени", null, tracesChart));
  },
  update(d, ctx) {
    data = d;
    const p = d.params;
    timesInput.value = [0, 0.25, 0.5, 0.75, 1].map((f) => parseFloat((f * p.T).toPrecision(6))).join(", ");
    pointsInput.value = [0.25, 0.5, 0.75].map((f) => parseFloat((f * p.l).toPrecision(6))).join(", ");
    load(ctx);
  },
};
