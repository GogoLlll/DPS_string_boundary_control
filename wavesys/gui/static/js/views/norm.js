// Вкладка «Норма и энергия»: насколько состояние струны далеко от покоя.
import { h, clear, card } from "../dom.js";
import { colors, hexToRgba, layout, line, plot } from "../theme.js";

let phiChart, energyBox;

export default {
  id: "norm",
  title: "Норма и энергия",
  mount(root) {
    phiChart = h("div", { class: "chart short" });
    energyBox = h("div", {});
    root.append(
      card("Норма состояния Φ(t)",
        "Φ(t) = ∫ y² dx + (1/a²) min по q ∫ (Q − q)² dx, где Q — первообразная скорости по x. Равна нулю тогда и только тогда, когда струна покоится. Конечна для любых данных, в том числе при скачках.",
        phiChart),
      card("Энергия E(t)", "E(t) = ½ ∫ (y_x² + y_t² / a²) dx. Без управления сохраняется; управление может её увеличивать и уменьшать.", energyBox));
  },
  update(data) {
    const c = colors(), n = data.norm, p = data.params;
    plot(phiChart, [line(c, c.s1, { x: n.t, y: n.phi, fill: "tozeroy", fillcolor: hexToRgba(c.s1, 0.1), hovertemplate: "Φ = %{y:.4g}<extra></extra>" })],
      layout(c, { x: "t", y: "Φ", xaxis: { range: [0, p.T] }, yaxis: { rangemode: "tozero" } }));
    clear(energyBox);
    if (n.energy) {
      const el = h("div", { class: "chart short" });
      energyBox.append(el);
      plot(el, [line(c, c.s1, { x: n.t, y: n.energy, fill: "tozeroy", fillcolor: hexToRgba(c.s1, 0.1), hovertemplate: "E = %{y:.4g}<extra></extra>" })],
        layout(c, { x: "t", y: "E", xaxis: { range: [0, p.T] }, yaxis: { rangemode: "tozero" } }));
    } else {
      energyBox.append(h("p", { class: "empty" },
        "В решении есть скачки, поэтому энергия бесконечна и не показывается. Близость к покою в этом случае измеряет норма Φ(t)."));
    }
  },
};
