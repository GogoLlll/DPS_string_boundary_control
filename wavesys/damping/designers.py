import numpy as np

from .. import design
from . import linking
from .base_profile import BaseProfile
from .controls import ContinuousControl, EqualShareControl

AFTER = ("наблюдать после T, в долях t₀", "0.5")
R_OPTION = ("параметр r (пусто — наилучший)", "")

def _number(text, name, default=None):
    if text is None or str(text).strip() == "":
        if default is None:
            raise ValueError(f"Не задан параметр «{name}»")
        return default
    try:
        value = float(str(text).replace(",", "."))
    except ValueError:
        raise ValueError(f"Параметр «{name}» должен быть числом, получено: {text!r}") from None
    if not np.isfinite(value):
        raise ValueError(f"Параметр «{name}» должен быть конечным числом")
    return value

def regime(profile):
    t0 = profile.params.t0
    if profile.kT == 2 and profile.sT == 0.0:
        return "критический", "T = 2·t₀"
    if profile.supercritical:
        return "сверхкритический", f"T = {profile.T / t0:.4g}·t₀ > 2·t₀"
    return "докритический", f"T = {profile.T / t0:.4g}·t₀ < 2·t₀"

def attach_report(control, params, data, mode, notes=()):
    T = control.T
    profile = BaseProfile(params, data, T)
    name, relation = regime(profile)
    phi, _ = linking.residual(control, params, data, T)
    norm0 = profile.initial_norm()
    exact = phi <= 1e-10 * max(norm0, 1e-300)
    sigma, G, n = profile.sample(400)
    y00 = float(np.asarray(data.y0(np.array([0.0]))).ravel()[0])
    y0l = float(np.asarray(data.y0(np.array([params.l]))).ravel()[0])
    notes = list(notes)
    if not profile.supercritical:
        notes.append("Время меньше критического 2·t₀: управление определено однозначно, свободы выбора нет."
                     if exact else
                     "Время меньше критического 2·t₀, а начальные данные общего вида: точное успокоение невозможно. "
                     "Построено управление с наименьшим остатком.")
    control.damping = {
        "mode": mode, "regime": name, "relation": relation, "target_time": T,
        "r": control.r, "r_optimal": profile.r_optimal(), "r_continuous": y00,
        "cost": linking.cost(control, params, data, T), "cost_min": profile.cost_min(),
        "residual": phi, "residual_min": profile.residual_min(), "initial_norm": norm0, "exact": bool(exact),
        "minimal_time": profile.minimal_time(), "y0_at_right_end": y0l,
        "sigma": sigma, "G": G, "count": n, "target": G + control.r / 2,
        "rho": profile.rho, "K": profile.K, "notes": notes,
    }
    return control

class MinimumCost(design.ControlDesigner):
    name = "Успокоение: наименьшая стоимость"
    description = ("Управление с наименьшим ∫u² dt, приводящее струну в покой к моменту T. "
                   "При T < 2·t₀ — управление с наименьшим остатком.")
    options = {"r": R_OPTION, "after": AFTER}

    def build(self, params, data, horizon, r="", **_):
        profile = BaseProfile(params, data, horizon)
        value = _number(r, "r", default=profile.r_optimal())
        notes = []
        if str(r).strip() != "" and abs(value - profile.r_optimal()) > 1e-12:
            notes.append("Параметр r задан вручную и отличается от наилучшего."
                         + ("" if profile.supercritical else " При T < 2·t₀ это увеличивает остаток."))
        control = EqualShareControl(params, data, horizon, value)
        control.label = "наименьшая стоимость: (G + r/2) / n"
        return attach_report(control, params, data, "наименьшая стоимость", notes)

class Continuous(design.ControlDesigner):
    name = "Успокоение: непрерывное управление"
    description = ("Непрерывное управление, приводящее струну в покой к моменту T ≥ 2·t₀. "
                   "Чем уже переходный слой, тем ближе стоимость к наименьшей и тем круче переходы.")
    options = {"eps": ("ширина переходного слоя, в долях t₀", "0.1"), "r": R_OPTION, "after": AFTER}

    def build(self, params, data, horizon, eps="0.1", r="", **_):
        profile = BaseProfile(params, data, horizon)
        if not profile.supercritical:
            raise ValueError("При T < 2·t₀ управление определено однозначно и выбрать непрерывное нельзя — "
                             "используйте режим «наименьшая стоимость».")
        notes = []
        y0l = float(np.asarray(data.y0(np.array([params.l]))).ravel()[0])
        y00 = float(np.asarray(data.y0(np.array([0.0]))).ravel()[0])
        if abs(y0l) > 1e-9 * max(1.0, abs(y00)):
            notes.append("y₀(l) ≠ 0: базовый профиль разрывен при σ = t₀, поэтому управление останется разрывным "
                         "в моменты t₀ + 2m·t₀. Непрерывность возможна только при y₀(l) = 0.")
        if profile.kT == 2 and profile.sT == 0.0:
            control = EqualShareControl(params, data, horizon, y00)
            control.label = "согласованное управление: G + y₀(0)/2"
            notes.append("При T = 2·t₀ непрерывное управление единственно: r = y₀(0); переходный слой не нужен.")
            return attach_report(control, params, data, "непрерывное", notes)
        width = _number(eps, "ширина переходного слоя", default=0.1)
        if not width > 0:
            raise ValueError("Ширина переходного слоя должна быть положительной")
        value = _number(r, "r", default=profile.r_optimal())
        control = ContinuousControl(params, data, horizon, value, width * params.t0)
        control.label = "непрерывное: равные доли с переходными слоями"
        if control.eps < width * params.t0 * (1 - 1e-12):
            notes.append(f"Слой сужен до {control.eps / params.t0:.4g}·t₀, чтобы слои у соседних точек не перекрывались.")
        return attach_report(control, params, data, "непрерывное", notes)

class TwoPassesThenRest(design.ControlDesigner):
    name = "Успокоение за 2·t₀, затем покой"
    description = "Управление лекции 6 на (0, 2·t₀), после чего левый конец неподвижен. Время T должно быть не меньше 2·t₀."
    options = {"variant": ("вариант", "оптимальное h", ["оптимальное h", "согласованное (непрерывное)"]), "after": AFTER}

    def build(self, params, data, horizon, variant="оптимальное h", **_):
        if horizon < 2 * params.t0 * (1 - 1e-12):
            raise ValueError("Этот режим требует T ≥ 2·t₀")
        T = 2 * params.t0
        profile = BaseProfile(params, data, T)
        if str(variant).startswith("согласованное"):
            r = float(np.asarray(data.y0(np.array([0.0]))).ravel()[0])
            label = "согласованное управление на (0, 2·t₀): G + y₀(0)/2"
        else:
            r = profile.r_optimal()
            label = "оптимальное управление лекции 6 на (0, 2·t₀): h"
        control = EqualShareControl(params, data, T, r)
        control.label = label
        notes = ["Покой наступает в момент 2·t₀ и сохраняется; стоимость не зависит от T."]
        return attach_report(control, params, data, "за 2·t₀, затем покой", notes)

for _designer in (MinimumCost(), Continuous(), TwoPassesThenRest()):
    design.register(_designer)
