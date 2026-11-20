// Цвета берутся из стилей страницы, чтобы графики следовали светлой и тёмной теме.
export function colors() {
  const cs = getComputedStyle(document.documentElement);
  const v = (name) => cs.getPropertyValue(name).trim();
  return {
    surface: v("--surface"), ink: v("--ink"), ink2: v("--ink-2"), muted: v("--muted"),
    grid: v("--grid"), axis: v("--axis"),
    s1: v("--series-1"), s2: v("--series-2"), reference: v("--reference"),
    ordinal: [1, 2, 3, 4, 5].map((i) => v(`--ord-${i}`)),
    diverging: ["--div-n4", "--div-n3", "--div-n2", "--div-n1", "--div-0", "--div-p1", "--div-p2", "--div-p3", "--div-p4"].map(v),
  };
}

// Расходящаяся шкала: синий — отрицательные значения, серый — ноль, красный — положительные.
export function divergingScale(c) {
  const n = c.diverging.length - 1;
  return c.diverging.map((col, i) => [i / n, col]);
}

// n оттенков одного цвета для упорядоченных кривых (n <= 5)
export function ordinalColors(c, n) {
  if (n <= 1) return [c.ordinal[2]];
  return Array.from({ length: n }, (_, i) => c.ordinal[Math.round((i * 4) / (n - 1))]);
}

export function hexToRgba(hex, alpha) {
  const m = /^#?([0-9a-f]{6})$/i.exec(hex);
  if (!m) return hex;
  const n = parseInt(m[1], 16);
  return `rgba(${n >> 16}, ${(n >> 8) & 255}, ${n & 255}, ${alpha})`;
}

function axis(c, title, extra = {}) {
  return {
    title: { text: title, font: { color: c.ink2, size: 12 }, standoff: 8 },
    gridcolor: c.grid, gridwidth: 1, linecolor: c.axis, zerolinecolor: c.axis, zerolinewidth: 1,
    ticks: "outside", tickcolor: c.axis, ticklen: 4, tickfont: { color: c.muted, size: 11 },
    automargin: true, ...extra,
  };
}

export function layout(c, { x = "", y = "", legend = false, xaxis = {}, yaxis = {}, ...rest } = {}) {
  return {
    paper_bgcolor: c.surface, plot_bgcolor: c.surface,
    font: { family: 'system-ui, -apple-system, "Segoe UI", sans-serif', color: c.ink2, size: 12 },
    margin: { l: 56, r: 16, t: legend ? 34 : 12, b: 46 },
    hovermode: "x unified",
    hoverlabel: { bgcolor: c.surface, bordercolor: c.axis, font: { color: c.ink, size: 12 } },
    showlegend: legend,
    legend: { orientation: "h", x: 0, y: 1.02, yanchor: "bottom", font: { color: c.ink2, size: 12 }, bgcolor: "rgba(0,0,0,0)" },
    xaxis: axis(c, x, xaxis), yaxis: axis(c, y, yaxis),
    ...rest,
  };
}

const CONFIG = {
  displaylogo: false, responsive: true,
  modeBarButtonsToRemove: ["lasso2d", "select2d", "autoScale2d", "toggleSpikelines", "hoverClosestCartesian", "hoverCompareCartesian"],
};

export function plot(el, traces, lay) {
  return window.Plotly.react(el, traces, lay, CONFIG);
}

export function line(c, color, extra = {}) {
  return { type: "scatter", mode: "lines", line: { color, width: 2, shape: "linear" }, connectgaps: false, ...extra };
}
