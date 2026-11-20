// Вкладка «Таблица»: значения y(t, x) на сетке и выгрузка в CSV.
import { getJSON } from "../api.js";
import { h, clear, card } from "../dom.js";
import { fmt } from "../format.js";

let rowsSel, colsSel, wrap, link, message;

function options(values, chosen) {
  return values.map((v) => h("option", { value: String(v), selected: v === chosen }, String(v)));
}

async function load() {
  message.textContent = "";
  const params = { rows: rowsSel.value, cols: colsSel.value };
  link.href = `/api/table.csv?${new URLSearchParams(params)}`;
  try {
    const t = await getJSON("/api/table", params);
    const head = h("tr", {}, h("th", { class: "num" }, "t \\ x"), ...t.x.map((x) => h("th", { class: "num" }, fmt(x, 5))));
    const body = t.t.map((time, i) => h("tr", {}, h("th", { class: "num" }, fmt(time, 5)), ...t.Y[i].map((v) => h("td", { class: "num" }, fmt(v, 5)))));
    clear(wrap).append(h("table", { class: "data" }, h("thead", {}, head), h("tbody", {}, ...body)));
  } catch (err) {
    message.textContent = err.message;
  }
}

export default {
  id: "table",
  title: "Таблица",
  mount(root) {
    rowsSel = h("select", { onchange: load }, ...options([5, 9, 13, 21, 41, 81], 13));
    colsSel = h("select", { onchange: load }, ...options([5, 6, 9, 11, 21, 41], 11));
    link = h("a", { href: "/api/table.csv", download: "values.csv" }, "Скачать CSV");
    wrap = h("div", { class: "scroll" });
    message = h("p", { class: "hint", role: "alert" });
    root.append(card("Значения y(t, x)", "Строки — моменты времени, столбцы — точки струны. Значения точные, а не взятые с сетки графиков.",
      h("div", { class: "controls" }, h("label", {}, "строк (по t)", rowsSel), h("label", {}, "столбцов (по x)", colsSel), link), message, wrap));
  },
  update() { load(); },
};
