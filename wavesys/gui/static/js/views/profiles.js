// Вкладка «Профили»: бегущие волны A(z) и B(z), из которых складывается решение.
import { h, card } from "../dom.js";
import { breakAt, profileJumps } from "../fronts.js";
import { colors, layout, line, plot } from "../theme.js";

let chart;

export default {
  id: "profiles",
  title: "Профили A, B",
  mount(root) {
    chart = h("div", { class: "chart tall" });
    const sub = h("span", { html: "Решение равно y(t, x) = ½ [A(t + x/a) + B(t + (l − x)/a)]. A — волна, бегущая влево, B — вправо. "
      + "Вертикальные линии — границы участков длины t₀; с участка на участок профили переходят по правилу "
      + "A<sub>k+1</sub> = −B<sub>k</sub>, &nbsp;B<sub>k+1</sub> = 2u<sub>k</sub> − A<sub>k</sub>." });
    root.append(card("Профили бегущих волн", sub, chart));
  },
  update(data) {
    const c = colors(), p = data.params, pr = data.profiles;
    const zMax = pr.z[pr.z.length - 1];
    const jA = [], jB = [];
    for (const src of data.sources) {
      const j = profileJumps(p, src, zMax);
      jA.push(...j.A); jB.push(...j.B);
    }
    const A = breakAt(pr.z, pr.A, jA), B = breakAt(pr.z, pr.B, jB);
    const count = Math.floor(zMax / p.t0 + 1e-9);
    const shapes = [], tickvals = [], ticktext = [];
    for (let k = 0; k <= count; k++) {
      shapes.push({ type: "line", xref: "x", yref: "paper", x0: k * p.t0, x1: k * p.t0, y0: 0, y1: 1, line: { color: c.grid, width: 1 }, layer: "below" });
      tickvals.push(k * p.t0);
      ticktext.push(k === 0 ? "0" : k === 1 ? "t₀" : `${k} t₀`);
    }
    const sparse = count > 14;
    plot(chart, [
      line(c, c.s1, { x: A.x, y: A.y, name: "A(z) — волна влево", hovertemplate: "%{y:.4g}" }),
      line(c, c.s2, { x: B.x, y: B.y, name: "B(z) — волна вправо", hovertemplate: "%{y:.4g}" }),
    ], layout(c, { x: "аргумент профиля z", y: "", legend: true, shapes: sparse ? [] : shapes,
      xaxis: sparse ? { range: [0, zMax] } : { range: [0, zMax], tickvals, ticktext, showgrid: false } }));
  },
};
