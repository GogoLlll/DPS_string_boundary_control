// Построение элементов страницы. Текст всегда вставляется как текст, а не как разметка.
export function h(tag, attrs = {}, ...children) {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v === null || v === undefined || v === false) continue;
    if (k === "class") el.className = v;
    else if (k.startsWith("on") && typeof v === "function") el.addEventListener(k.slice(2), v);
    else if (k === "html") el.innerHTML = v;            // только для разметки, написанной в коде (значки)
    else el.setAttribute(k, v === true ? "" : v);
  }
  for (const c of children.flat()) {
    if (c === null || c === undefined || c === false) continue;
    el.append(c instanceof Node ? c : document.createTextNode(String(c)));
  }
  return el;
}

export function clear(el) {
  while (el.firstChild) el.removeChild(el.firstChild);
  return el;
}

const ICONS = {
  good: '<svg width="14" height="14" viewBox="0 0 14 14" aria-hidden="true"><circle cx="7" cy="7" r="7" fill="var(--good)"/><path d="M3.8 7.3l2.2 2.2 4.2-4.6" fill="none" stroke="#fff" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/></svg>',
  bad: '<svg width="14" height="14" viewBox="0 0 14 14" aria-hidden="true"><circle cx="7" cy="7" r="7" fill="var(--critical)"/><path d="M4.4 4.4l5.2 5.2M9.6 4.4l-5.2 5.2" stroke="#fff" stroke-width="1.8" stroke-linecap="round"/></svg>',
  info: '<svg width="14" height="14" viewBox="0 0 14 14" aria-hidden="true"><circle cx="7" cy="7" r="6.2" fill="none" stroke="var(--muted)" stroke-width="1.6"/><path d="M4.2 7h5.6" stroke="var(--muted)" stroke-width="1.6" stroke-linecap="round"/></svg>',
};

// Состояние всегда показывается значком и словом, а не одним цветом.
export function stateBadge(kind, text) {
  return h("span", { class: "state" }, h("span", { html: ICONS[kind] }), text);
}

export function card(title, subtitle, ...body) {
  return h("section", { class: "card" }, h("h2", {}, title), subtitle ? h("p", { class: "sub" }, subtitle) : null, ...body);
}
