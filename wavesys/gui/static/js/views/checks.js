// Вкладка «Проверки»: таблица всех проверок решения.
import { h, clear, card, stateBadge } from "../dom.js";
import { fmt, sci } from "../format.js";

let body, summary;

function value(v) {
  if (v === null) return "";
  return Math.abs(v) !== 0 && (Math.abs(v) < 1e-3 || Math.abs(v) >= 1e4) ? sci(v) : fmt(v);
}

export default {
  id: "checks",
  title: "Проверки",
  mount(root) {
    summary = h("p", { class: "sub" });
    body = h("tbody");
    root.append(card("Проверки решения", null, summary,
      h("div", { class: "scroll" }, h("table", { class: "data" },
        h("thead", {}, h("tr", {}, h("th", {}, "Что проверяется"), h("th", { class: "num" }, "Получено"), h("th", { class: "num" }, "Допуск"), h("th", {}, "Итог"), h("th", {}, "Пояснение"))),
        body))));
  },
  update(data) {
    const s = data.stats;
    summary.textContent = `Выполнено ${s.checks_passed} из ${s.checks_total} проверок с допуском. Строки «справочно» ничего не утверждают: они показывают, насколько гладким будет решение.`;
    clear(body);
    let group = null;
    for (const c of data.checks) {
      if (c.group !== group) {
        group = c.group;
        body.append(h("tr", { class: "group" }, h("td", { colspan: "5" }, group)));
      }
      const badge = c.ok === true ? stateBadge("good", "выполнено") : c.ok === false ? stateBadge("bad", "не выполнено") : stateBadge("info", "справочно");
      body.append(h("tr", {}, h("td", {}, c.name), h("td", { class: "num" }, value(c.value)),
        h("td", { class: "num" }, c.limit === null ? "" : `≤ ${value(c.limit)}`), h("td", {}, badge), h("td", {}, c.note || "")));
    }
  },
};
