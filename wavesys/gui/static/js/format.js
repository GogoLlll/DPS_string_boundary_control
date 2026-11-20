// Запись чисел. Очень малые и очень большие — в виде a·10^n.
const SUP = { "-": "⁻", 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹" };

export function sci(v, digits = 2) {
  if (v === null || v === undefined || Number.isNaN(v)) return "—";
  if (v === 0) return "0";
  const [m, e] = v.toExponential(digits).split("e");
  const exp = String(parseInt(e, 10)).split("").map((c) => SUP[c]).join("");
  return `${m}·10${exp}`.replace("-", "−");
}

export function fmt(v, digits = 4) {
  if (v === null || v === undefined || Number.isNaN(v)) return "—";
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a < 1e-4 || a >= 1e6) return sci(v, Math.max(1, digits - 2));
  return String(parseFloat(v.toPrecision(digits))).replace("-", "−");
}

// Время в долях пробега: «1.5 t₀»
export function inT0(t, t0) {
  return `${fmt(t / t0, 4)} t₀`;
}
