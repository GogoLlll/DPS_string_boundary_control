// Вкладка «Успокоение»: достигнута ли цель, сколько это стоит и из чего построено управление.
import { h, clear, card, stateBadge } from "../dom.js";
import { fmt, inT0 } from "../format.js";
import { breakAt } from "../fronts.js";
import { colors, layout, line, plot } from "../theme.js";

let box;

function tile(label, value, foot, small = false) {
  return h("div", { class: "tile" }, h("div", { class: "label" }, label), h("div", { class: small ? "value small" : "value" }, value), h("div", { class: "foot" }, foot));
}

function row(name, value, note) {
  return h("tr", {}, h("td", {}, name), h("td", { class: "num" }, value), h("td", {}, note || ""));
}

export default {
  id: "damping",
  title: "Успокоение",
  mount(root) {
    box = h("div", {});
    root.append(box);
  },
  update(data) {
    clear(box);
    const d = data.damping;
    if (!d) {
      box.append(card("Задача успокоения", null, h("p", { class: "empty" },
        "Сейчас управление задано формулой. Чтобы программа сама построила управление, приводящее струну в покой к моменту T, "
        + "выберите в поле «Источник» один из режимов «Успокоение…» и нажмите «Рассчитать».")));
      return;
    }
    const c = colors(), p = data.params, t0 = p.t0, T = d.target_time;
    const best = d.cost_min === null ? null : d.cost_min;

    box.append(h("div", { class: "tiles" },
      tile("Цель", d.exact ? stateBadge("good", "покой достигнут") : stateBadge("bad", "покой недостижим"),
        `к моменту T = ${fmt(T)} (${inT0(T, t0)})`, true),
      tile("Норма состояния в момент T", fmt(d.state_norm_at_T), d.exact ? "Φ(T); ноль — покой" : "Φ(T) — наименьший возможный остаток"),
      tile("Стоимость управления", fmt(d.cost), "J = ∫ u² dt по (0, T)"),
      tile("Наименьшая возможная стоимость", best === null ? "—" : fmt(best),
        best === null ? "при T < 2·t₀ управление единственно" : `J_min(T); превышение ${fmt(100 * (d.cost / best - 1), 3)} %`),
      tile("Режим", d.regime, d.relation, true),
      tile("Минимальное время для этих данных", `${fmt(d.minimal_time / t0, 4)} t₀`, "раньше покой недостижим"),
    ));

    // правая часть условия связи: то, что управление должно набрать в сумме по повторам
    const chart = h("div", { class: "chart" });
    const legend = h("div", { class: "legend" });
    box.append(card("Условие связи: что должно набрать управление",
      "Сумма значений управления в моменты σ, σ + 2t₀, σ + 4t₀, … из (0, T) должна равняться G(σ) + r/2 — линия на графике. "
      + "Число n(σ) таких моментов указано над графиком; где n = 0, управление до этой части не дотягивается.",
      legend, chart));
    const cuts = [t0];
    if (d.rho > 0) cuts.push(d.rho);
    const target = breakAt(d.sigma, d.target, cuts);
    const shapes = [], annotations = [];
    const marks = [0, ...(d.rho > 0 && d.rho < 2 * t0 ? [d.rho] : []), 2 * t0];
    if (d.rho > 0) shapes.push({ type: "line", xref: "x", yref: "paper", x0: d.rho, x1: d.rho, y0: 0, y1: 1, line: { color: c.axis, width: 1 } });
    shapes.push({ type: "line", xref: "x", yref: "paper", x0: t0, x1: t0, y0: 0, y1: 1, line: { color: c.grid, width: 1 } });
    for (let i = 0; i + 1 < marks.length; i++) {
      const mid = 0.5 * (marks[i] + marks[i + 1]);
      const idx = d.sigma.findIndex((v) => v >= mid);
      annotations.push({ x: mid, y: 1.02, xref: "x", yref: "paper", yanchor: "bottom", showarrow: false,
        text: `n = ${d.count[Math.max(idx, 0)]}`, font: { color: c.ink2, size: 12 } });
    }
    plot(chart, [line(c, c.s1, { x: target.x, y: target.y, hovertemplate: "G + r/2 = %{y:.4g}<extra></extra>" })],
      layout(c, { x: "σ — момент прихода возмущения к левому концу", y: "G(σ) + r/2", shapes, annotations, margin: { l: 56, r: 16, t: 30, b: 46 },
        xaxis: { range: [0, 2 * t0], tickvals: [0, t0, 2 * t0], ticktext: ["0", "t₀", "2 t₀"] } }));

    const rows = [
      row("Способ построения", d.mode, data.texts.u),
      row("Момент успокоения T", fmt(T), `${inT0(T, t0)}; T = 2K·t₀ + ρ, K = ${d.K}, ρ = ${fmt(d.rho)}`),
      row("Параметр r", fmt(d.r, 6), "сдвиг базового профиля: управление набирает G + r/2"),
      row(d.regime === "докритический" ? "r с наименьшим остатком" : "r с наименьшей стоимостью", fmt(d.r_optimal, 6),
        d.regime === "докритический" ? "−2 × среднее G на (T, 2t₀)" : "r*(T) = −2M/W; при T = 2K·t₀ равно −F̄/a"),
      row("r для согласования в начале", fmt(d.r_continuous, 6), "r = y₀(0): управление начинается со значения y₀(0)"),
      row("Стоимость J", fmt(d.cost, 6), "∫ u² dt по (0, T)"),
      row("Наименьшая стоимость J_min(T)", best === null ? "—" : fmt(best, 6), best === null ? "определена при T ≥ 2·t₀" : "S − M²/W"),
      row("Норма начального состояния Φ(0)", fmt(d.initial_norm, 6), "расстояние до покоя в начале"),
      row("Остаток по управлению", fmt(d.residual, 6), "Φ(T), вычисленное только из управления (условие связи)"),
      row("Остаток по состоянию", fmt(d.state_norm_at_T, 6), "Φ(T), вычисленное по движению струны"),
      row("Наименьший возможный остаток Φ_min(T)", fmt(d.residual_min, 6), "ноль при T ≥ 2·t₀"),
      row("Норма в конце наблюдения", fmt(d.state_norm_at_end, 6), `t = ${fmt(p.T)}; после T управление выключено`),
      row("Минимальное время T_min", fmt(d.minimal_time, 5), `${inT0(d.minimal_time, t0)}; оценка по сетке`),
    ];
    box.append(card("Числа", null, h("div", { class: "scroll" }, h("table", { class: "data" },
      h("thead", {}, h("tr", {}, h("th", {}, "Величина"), h("th", { class: "num" }, "Значение"), h("th", {}, "Пояснение"))),
      h("tbody", {}, ...rows)))));

    const notes = h("ul", { class: "notes" });
    for (const n of d.notes) notes.append(h("li", {}, n));
    if (d.regime !== "докритический" && Math.abs(d.r - d.r_continuous) > 1e-9 && d.mode !== "непрерывное") {
      notes.append(h("li", {}, "Управление начинается не со значения y₀(0), поэтому по струне бежит скачок; к моменту T все скачки гасятся. "
        + "Режим «непрерывное управление» избавляет от скачков ценой большей стоимости."));
    }
    notes.append(h("li", {}, "Формулы и их вывод — в reshenie.pdf: условие связи (XIX), наименьшая стоимость (XXIII)–(XXIV), остаток (XXXV), непрерывное управление (XLII)."));
    box.append(card("Замечания", null, notes));
  },
};
