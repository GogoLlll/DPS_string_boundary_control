// Сборка страницы: форма, вкладки, запуск расчёта.
import { getJSON, postJSON } from "./api.js";
import { h } from "./dom.js";
import { fmt } from "./format.js";
import { initForm } from "./form.js";
import overview from "./views/overview.js";
import motion from "./views/motion.js";
import sections from "./views/sections.js";
import control from "./views/control.js";
import profiles from "./views/profiles.js";
import norm from "./views/norm.js";
import checks from "./views/checks.js";
import table from "./views/table.js";
import point from "./views/point.js";
import damping from "./views/damping.js";

// Чтобы добавить вкладку, достаточно написать модуль в views/ и вписать его сюда.
const VIEWS = [overview, motion, damping, sections, control, profiles, norm, checks, table, point];

const state = { data: null, active: VIEWS[0].id, form: null, dirty: new Set() };
const tabsEl = document.getElementById("tabs");
const viewsEl = document.getElementById("views");
const statusEl = document.getElementById("status");
const panels = new Map();

function setStatus(text, isError = false) {
  statusEl.textContent = text;
  statusEl.classList.toggle("error", isError);
}

const ctx = {
  data: () => state.data,
  goto(id, args) {
    activate(id);
    const view = VIEWS.find((v) => v.id === id);
    if (view.receive) view.receive(args, ctx);
  },
  setStatus,
};

function refresh(view) {
  if (!state.data || !state.dirty.has(view.id)) return;
  state.dirty.delete(view.id);
  try {
    view.update(state.data, ctx);
  } catch (err) {
    console.error(err);
    setStatus(`Не удалось показать вкладку «${view.title}»: ${err.message}`, true);
  }
}

function activate(id) {
  const prev = VIEWS.find((v) => v.id === state.active);
  if (prev && prev.id !== id && prev.hide) prev.hide();
  state.active = id;
  const q = new URLSearchParams(window.location.hash.replace(/^#/, ""));
  q.set("tab", id);
  history.replaceState(null, "", `#${q}`);
  for (const view of VIEWS) {
    const on = view.id === id;
    panels.get(view.id).hidden = !on;
    const tab = document.getElementById(`tab-${view.id}`);
    tab.setAttribute("aria-selected", on ? "true" : "false");
    tab.tabIndex = on ? 0 : -1;
  }
  const view = VIEWS.find((v) => v.id === id);
  refresh(view);
  if (view.show) view.show(ctx);
  window.dispatchEvent(new Event("resize"));           // графики подстраиваются под ширину вкладки
}

function buildTabs() {
  for (const view of VIEWS) {
    const tab = h("button", { class: "tab", id: `tab-${view.id}`, role: "tab", type: "button", "aria-controls": `view-${view.id}`,
      onclick: () => activate(view.id) }, view.title);
    tabsEl.append(tab);
    const panel = h("div", { class: "view", id: `view-${view.id}`, role: "tabpanel", "aria-labelledby": `tab-${view.id}`, hidden: true });
    viewsEl.append(panel);
    panels.set(view.id, panel);
    view.mount(panel, ctx);
  }
  tabsEl.addEventListener("keydown", (ev) => {
    if (ev.key !== "ArrowRight" && ev.key !== "ArrowLeft") return;
    const i = VIEWS.findIndex((v) => v.id === state.active);
    const next = VIEWS[(i + (ev.key === "ArrowRight" ? 1 : VIEWS.length - 1)) % VIEWS.length];
    activate(next.id);
    document.getElementById(`tab-${next.id}`).focus();
  });
}

async function run(payload) {
  state.form.setBusy(true);
  viewsEl.classList.add("stale");                       // прежние графики остаются на месте, пока идёт расчёт
  setStatus("Идёт расчёт…");
  const started = performance.now();
  try {
    state.data = await postJSON("/api/solve", payload);
    for (const view of VIEWS) state.dirty.add(view.id);
    const s = state.data.stats;
    const failed = s.checks_total - s.checks_passed;
    setStatus(`Готово за ${fmt((performance.now() - started) / 1000, 2)} с. Проверок выполнено: ${s.checks_passed} из ${s.checks_total}`
      + (failed ? ` — ${failed} не выполнено, см. вкладку «Проверки».` : "."), failed > 0);
    activate(state.active);
  } catch (err) {
    setStatus(err.message, true);
  } finally {
    viewsEl.classList.remove("stale");
    state.form.setBusy(false);
  }
}

// Адрес вида  #tab=motion&preset=2  открывает нужную вкладку и пример.
function fromHash() {
  const q = new URLSearchParams(window.location.hash.replace(/^#/, ""));
  return { tab: q.get("tab"), preset: q.get("preset") };
}

async function init() {
  const start = fromHash();
  buildTabs();
  activate(VIEWS.some((v) => v.id === start.tab) ? start.tab : state.active);
  try {
    const meta = await getJSON("/api/meta");
    state.form = initForm(meta, run, start.preset);
    await run(state.form.collect());
  } catch (err) {
    setStatus(err.message, true);
  }
}

// Ошибки самой страницы показываются пользователю, а не теряются в консоли браузера.
window.addEventListener("error", (ev) => setStatus(`Сбой интерфейса: ${ev.message}`, true));
window.addEventListener("unhandledrejection", (ev) => setStatus(`Сбой интерфейса: ${ev.reason && ev.reason.message ? ev.reason.message : ev.reason}`, true));

init();
