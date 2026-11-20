// Форма с параметрами задачи: примеры, поля, сбор данных для расчёта.
import { h, clear } from "./dom.js";
import { fmt } from "./format.js";

const $ = (id) => document.getElementById(id);
const FIELDS = { l: "f-l", a: "f-a", T: "f-T", y0: "f-y0", y1: "f-y1", u: "f-u" };

export function initForm(meta, onSubmit, startPreset) {
  const preset = $("preset");
  for (const name of Object.keys(meta.presets)) preset.append(h("option", { value: name }, name));
  preset.append(h("option", { value: "" }, "— свои данные —"));

  const source = $("f-source");
  for (const s of meta.sources) source.append(h("option", { value: s.name }, s.name));

  const help = $("help-list");
  for (const [what, means] of meta.functions) {
    help.append(h("dt", {}, what));
    if (means) help.append(h("dd", {}, means));
  }

  function applyPreset(name) {
    const p = meta.presets[name];
    if (!p) { $("preset-note").textContent = ""; return; }
    for (const [key, id] of Object.entries(FIELDS)) $(id).value = p[key];
    source.value = p.source && meta.sources.some((v) => v.name === p.source) ? p.source : meta.sources[0].name;
    $("preset-note").textContent = p.note || "";
    updateSource();
    for (const [key, value] of Object.entries(p.options || {})) {
      const el = document.querySelector(`[data-option="${key}"]`);
      if (el) el.value = value;
    }
    updateHint();
  }

  function updateHint() {
    const l = parseFloat($("f-l").value.replace(",", ".")), a = parseFloat($("f-a").value.replace(",", "."));
    const T = parseFloat($("f-T").value.replace(",", "."));
    $("t0-hint").textContent = l > 0 && a > 0
      ? `Время пробега t₀ = l / a = ${fmt(l / a)}` + (T > 0 ? `; T = ${fmt(T / (l / a))} t₀` : "")
      : "";
  }

  function updateSource() {
    const s = meta.sources.find((v) => v.name === source.value) || meta.sources[0];
    $("source-note").textContent = s.description || "";
    $("u-field").hidden = s !== meta.sources[0];
    const box = clear($("source-options"));
    for (const [key, opt] of Object.entries(s.options || {})) {
      const input = opt.choices
        ? h("select", { "data-option": key }, ...opt.choices.map((c) => h("option", { value: c, selected: c === opt.default }, c)))
        : h("input", { "data-option": key, value: opt.default });
      input.addEventListener("input", clearPreset);
      box.append(h("label", { class: "field" }, h("span", {}, opt.label), input));
    }
  }

  function clearPreset() { preset.value = ""; $("preset-note").textContent = ""; }

  preset.addEventListener("change", () => applyPreset(preset.value));
  source.addEventListener("change", () => { clearPreset(); updateSource(); });
  for (const id of ["f-l", "f-a", "f-T"]) $(id).addEventListener("input", updateHint);
  for (const id of Object.values(FIELDS)) $(id).addEventListener("input", clearPreset);

  $("problem-form").addEventListener("submit", (ev) => {
    ev.preventDefault();
    onSubmit(collect());
  });

  const names = Object.keys(meta.presets);
  const first = names[parseInt(startPreset, 10)] || (names.includes(startPreset) ? startPreset : meta.default);
  preset.value = first;
  applyPreset(first);
  return { collect, setBusy: (busy) => { $("run").disabled = busy; $("run").textContent = busy ? "Считаю…" : "Рассчитать"; } };
}

export function collect() {
  const payload = { source: $("f-source").value, nx: $("f-nx").value, nt: $("f-nt").value || "auto", fd_nodes: $("f-fd").value, options: {} };
  for (const [key, id] of Object.entries(FIELDS)) payload[key] = $(id).value;
  for (const el of document.querySelectorAll("[data-option]")) payload.options[el.dataset.option] = el.value;
  return payload;
}
