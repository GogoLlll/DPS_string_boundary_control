// Вкладка «Обзор»: главные числа и карта решения на плоскости (x, t).
import { h, clear, card, stateBadge } from "../dom.js";
import { fmt, inT0 } from "../format.js";
import { colors, divergingScale, layout, plot } from "../theme.js";

let tiles, chart, notes, legend, clickBound = false;

function tile(label, value, foot) {
  return h("div", { class: "tile" }, h("div", { class: "label" }, label), h("div", { class: "value" }, value), h("div", { class: "foot" }, foot));
}

export default {
  id: "overview",
  title: "Обзор",
  mount(root) {
    tiles = h("div", { class: "tiles" });
    chart = h("div", { class: "chart tall" });
    legend = h("div", { class: "legend" });
    notes = h("ul", { class: "notes" });
    root.append(
      tiles,
      card("Решение y(t, x)", "По горизонтали — точка струны x, по вертикали — время t. Синий — отклонение вниз, красный — вверх. Щелчок по карте открывает разбор точки.",
        legend, chart),
      card("Что получилось", null, notes));
  },
  update(data, ctx) {
    const c = colors(), s = data.stats, p = data.params;
    clear(tiles).append(
      tile("Время пробега", fmt(s.t0), "t₀ = l / a"),
      tile("Время наблюдения", `${fmt(s.T_over_t0)} t₀`, `до t = ${fmt(p.T)}`),
      tile("Наибольшее отклонение", fmt(s.amplitude), "max |y| на сетке"),
      tile("Норма состояния в начале", fmt(s.phi_start), "Φ(0)"),
      tile("Норма состояния в конце", fmt(s.phi_end), "Φ в конце наблюдения; ноль — покой"),
      tile("Проверки", `${s.checks_passed} из ${s.checks_total}`,
        s.checks_passed === s.checks_total ? stateBadge("good", "все выполнены") : stateBadge("bad", "есть невыполненные")),
    );

    const M = s.amplitude > 0 ? s.amplitude : 1;
    const traces = [{
      type: "heatmap", x: data.grid.x, y: data.grid.t, z: data.grid.Y, zmin: -M, zmax: M, zmid: 0,
      colorscale: divergingScale(c), zsmooth: false,
      colorbar: { title: { text: "y", font: { color: c.ink2 } }, thickness: 12, outlinewidth: 0, tickfont: { color: c.muted, size: 11 }, len: 0.9 },
      hovertemplate: "x = %{x:.4g}<br>t = %{y:.4g}<br>y = %{z:.4g}<extra></extra>",
    }];
    for (const f of data.fronts) {
      traces.push({ type: "scatter", mode: "lines", x: f.x, y: f.t, line: { color: c.ink, width: 1.5 }, hoverinfo: "skip", showlegend: false });
    }
    clear(legend);
    if (data.fronts.length) {
      const key = h("span", { class: "swatch" });
      key.style.borderTopColor = c.ink;
      legend.append(h("span", { class: "key" }, key, "линия разрыва решения (фронт)"));
    }
    const goal = data.damping && data.damping.target_time < p.T * (1 - 1e-9) ? data.damping.target_time : null;
    if (goal !== null) {
      const key = h("span", { class: "swatch" });
      key.style.borderTopColor = c.muted;
      legend.append(h("span", { class: "key" }, key, `момент успокоения T = ${fmt(goal)}`));
    }
    const shapes = goal === null ? [] : [{ type: "line", xref: "paper", yref: "y", x0: 0, x1: 1, y0: goal, y1: goal, line: { color: c.muted, width: 1 } }];
    plot(chart, traces, layout(c, { x: "x", y: "t", hovermode: "closest", shapes,
      xaxis: { range: [0, p.l], showgrid: false }, yaxis: { range: [0, p.T], showgrid: false }, margin: { l: 56, r: 8, t: 12, b: 46 } }));
    if (!clickBound) {
      chart.on("plotly_click", (ev) => {
        const pt = ev.points && ev.points[0];
        if (pt) ctx.goto("point", { t: pt.y, x: pt.x });
      });
      clickBound = true;
    }

    clear(notes);
    notes.append(h("li", {}, `Гладкость: ${s.verdict}.`));
    if (data.sources.length) {
      for (const src of data.sources) {
        notes.append(h("li", {}, `Скачок величины ${fmt(src.size)} — ${src.label}; фронт бежит со скоростью a и отражается от концов со сменой знака.`));
      }
    } else {
      notes.append(h("li", {}, "Скачков в решении нет."));
    }
    notes.append(h("li", {}, `Данные: y₀(x) = ${data.texts.y0}; y₁(x) = ${data.texts.y1}; u(t) = ${data.texts.u}.`));
    notes.append(h("li", {}, `Конец наблюдения: t = ${fmt(p.T)} = ${inT0(p.T, p.t0)}.`));
    if (data.damping) {
      const d = data.damping;
      notes.append(h("li", {}, d.exact
        ? `Цель успокоения достигнута: в момент T = ${fmt(d.target_time)} норма состояния Φ(T) = ${fmt(d.state_norm_at_T)}. Подробности — на вкладке «Успокоение».`
        : `Точное успокоение к моменту T = ${fmt(d.target_time)} невозможно; достигнут наименьший остаток Φ(T) = ${fmt(d.state_norm_at_T)}. Подробности — на вкладке «Успокоение».`));
    }
  },
};
