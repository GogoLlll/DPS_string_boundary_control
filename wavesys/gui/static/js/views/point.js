// Вкладка «Разбор точки»: из чего складывается значение y(t, x).
import { getJSON } from "../api.js";
import { h, clear, card } from "../dom.js";
import { fmt, inT0 } from "../format.js";

let tInput, xInput, out, data = null;

const REGIONS = {
  I: "I — значение определяется только начальными данными (формула Даламбера)",
  II: "II — дошло воздействие с левого конца",
  III: "III — дошло отражение от правого конца",
  IV: "IV — дошли оба воздействия",
};

function controlList(items) {
  if (!items.length) return "не входят";
  return items.map((m) => `u(${fmt(m.t)}) = ${fmt(m.u)}`).join(";  ");
}

async function load() {
  if (!data) return;
  clear(out);
  try {
    const r = await getJSON("/api/point", { t: tInput.value, x: xInput.value });
    const p = data.params;
    out.append(
      h("div", { class: "tiles" },
        h("div", { class: "tile" }, h("div", { class: "label" }, "Отклонение y(t, x)"), h("div", { class: "value" }, fmt(r.y, 6)), h("div", { class: "foot" }, "½ (A + B)")),
        h("div", { class: "tile" }, h("div", { class: "label" }, "Волна влево A"), h("div", { class: "value" }, fmt(r.A, 6)), h("div", { class: "foot" }, `участок A${r.kA}`)),
        h("div", { class: "tile" }, h("div", { class: "label" }, "Волна вправо B"), h("div", { class: "value" }, fmt(r.B, 6)), h("div", { class: "foot" }, `участок B${r.kB}`))),
      h("dl", { class: "kv" },
        h("dt", {}, "Точка"), h("dd", {}, `t = ${fmt(r.t, 6)} (${inT0(r.t, p.t0)}), x = ${fmt(r.x, 6)} (${fmt(r.x / p.l)} l)`),
        h("dt", {}, "Аргумент A"), h("dd", {}, `z = t + x/a = ${fmt(r.zA, 6)} = ${r.kA}·t₀ + ${fmt(r.sA, 6)}`),
        h("dt", {}, "Участок A"), h("dd", {}, `A${r.kA} строится из ${r.A_from === "A0" ? "A₀" : "−B₀"} и управления в моменты: ${controlList(r.control_A)}`),
        h("dt", {}, "Аргумент B"), h("dd", {}, `z = t + (l − x)/a = ${fmt(r.zB, 6)} = ${r.kB}·t₀ + ${fmt(r.sB, 6)}`),
        h("dt", {}, "Участок B"), h("dd", {}, `B${r.kB} строится из ${r.B_from === "B0" ? "B₀" : "−A₀"} и управления в моменты: ${controlList(r.control_B)}`),
        h("dt", {}, "Область первого слоя"), h("dd", {}, r.region ? REGIONS[r.region] : "точка лежит позже первого слоя (t ≥ t₀)"),
        h("dt", {}, "Функция Q(t, x)"), h("dd", {}, `${fmt(r.Q, 6)}  (a/2 · (A − B); её производная по x равна скорости)`)),
      h("p", { class: "hint" }, "Все моменты управления не позже t: значение в точке не зависит от будущего управления."));
  } catch (err) {
    out.append(h("p", { class: "hint", role: "alert" }, err.message));
  }
}

export default {
  id: "point",
  title: "Разбор точки",
  mount(root) {
    tInput = h("input", { type: "text", size: "10", "aria-label": "Момент t" });
    xInput = h("input", { type: "text", size: "10", "aria-label": "Точка x" });
    out = h("div", {});
    root.append(card("Из чего складывается значение в точке",
      "Значение y(t, x) — полусумма двух профилей. Здесь показано, на какие участки профилей попадает точка и какие значения управления в них входят.",
      h("form", { class: "controls", onsubmit: (ev) => { ev.preventDefault(); load(); } },
        h("label", {}, "t =", tInput), h("label", {}, "x =", xInput), h("button", { type: "submit" }, "Разобрать")),
      out));
  },
  update(d) {
    data = d;
    tInput.value = parseFloat((0.6 * d.params.T).toPrecision(6));
    xInput.value = parseFloat((0.4 * d.params.l).toPrecision(6));
    load();
  },
  receive(args) {
    tInput.value = parseFloat(args.t.toPrecision(6));
    xInput.value = parseFloat(args.x.toPrecision(6));
    load();
  },
};
