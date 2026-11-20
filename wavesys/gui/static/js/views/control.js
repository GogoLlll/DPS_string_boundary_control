// Вкладка «Управление»: смещение левого конца u(t).
import { h, clear, card } from "../dom.js";
import { fmt, inT0 } from "../format.js";
import { breakAt } from "../fronts.js";
import { colors, layout, line, plot } from "../theme.js";

let chart, list;

export default {
  id: "control",
  title: "Управление",
  mount(root) {
    chart = h("div", { class: "chart" });
    list = h("ul", { class: "notes" });
    root.append(card("Управление u(t) — смещение левого конца", null, chart), card("Особенности управления", null, list));
  },
  update(data) {
    const c = colors(), p = data.params, b = data.boundary;
    const jumps = data.sources.filter((s) => s.kind === "left" && s.time > 0);
    const cut = breakAt(b.t, b.u, jumps.map((s) => s.time));
    plot(chart, [line(c, c.s1, { x: cut.x, y: cut.y, hovertemplate: "u = %{y:.4g}<extra></extra>" })],
      layout(c, { x: "t", y: "u", xaxis: { range: [0, p.T] } }));
    clear(list);
    list.append(h("li", {}, `Формула: u(t) = ${data.texts.u}`));
    list.append(h("li", {}, `Значение в начале u(0) = ${fmt(b.u[0])}, в конце u(T) = ${fmt(b.u[b.u.length - 1])}.`));
    const corner = data.sources.find((s) => s.kind === "left" && s.time === 0);
    list.append(h("li", {}, corner
      ? `В начальный момент управление не согласовано с начальным отклонением: u(0) − y₀(0) = ${fmt(corner.size)}.`
      : "В начальный момент управление согласовано с начальным отклонением: u(0) = y₀(0)."));
    if (jumps.length) {
      for (const s of jumps) list.append(h("li", {}, `Скачок величины ${fmt(s.size)} при t = ${fmt(s.time)} (${inT0(s.time, p.t0)}).`));
    } else {
      list.append(h("li", {}, "Внутри интервала (0, T) скачков управления не обнаружено."));
    }
  },
};
