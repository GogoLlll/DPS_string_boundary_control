import json

import numpy as np
import pytest

from wavesys import damping, design, diagnostics  # noqa: F401
from wavesys.damping import designers
from wavesys.gui import api
from wavesys.presets import PRESETS
from wavesys.problem import Problem

MIN, CONT, TWO = "Успокоение: наименьшая стоимость", "Успокоение: непрерывное управление", "Успокоение за 2·t₀, затем покой"
BASE = dict(l="1", a="1", y0="0.8*(1 - x/l)", y1="0", u="0", nx=41, nt=41, fd_nodes=60)

def solve(T, source, **options):
    return api.solve(dict(BASE, T=str(T), source=source, options=options))

def failed(body):
    return [c["name"] for c in body["checks"] if c["ok"] is False]

def test_designers_are_registered_and_offered():
    assert {MIN, CONT, TWO} <= set(design.REGISTRY)
    sources = {s["name"]: s for s in api.meta()["sources"]}
    assert sources[TWO]["options"]["variant"]["choices"] == ["оптимальное h", "согласованное (непрерывное)"]
    assert sources[MIN]["options"]["after"]["default"] == "0.5"
    json.dumps(api.meta(), allow_nan=False)

def test_minimum_cost_at_critical_time_matches_lecture_6():
    r, body = solve(2, MIN, after="0.5")
    d = body["damping"]
    assert d["regime"] == "критический" and d["exact"] is True
    assert d["cost"] == pytest.approx(0.8 ** 2 / 6, rel=1e-9) and d["cost_min"] == pytest.approx(0.8 ** 2 / 6, rel=1e-9)
    assert d["r"] == pytest.approx(0.0, abs=1e-12) and d["r_continuous"] == pytest.approx(0.8)
    assert d["state_norm_at_T"] < 1e-25 and d["state_norm_at_end"] < 1e-25
    assert body["params"]["T"] == pytest.approx(2.5)                      # наблюдение продлено на 0.5 t0
    assert d["target_time"] == pytest.approx(2.0)
    assert failed(body) == []
    json.dumps(body, allow_nan=False)

def test_observation_extension_is_optional_and_validated():
    _, body = solve(2, MIN, after="")
    assert body["params"]["T"] == pytest.approx(2.0)
    with pytest.raises(api.RequestError):
        solve(2, MIN, after="-1")
    with pytest.raises(api.RequestError):
        solve(2, MIN, after="abc")

def test_manual_r_is_used_and_reported():
    _, body = solve(2, MIN, r="0.8", after="0")
    d = body["damping"]
    assert d["r"] == pytest.approx(0.8) and d["exact"] is True
    assert d["cost"] == pytest.approx(4 * 0.8 ** 2 / 6, rel=1e-9)          # согласованное управление вчетверо дороже
    assert any("вручную" in n for n in d["notes"])
    assert body["sources"] == []                                           # при r = y0(0) скачков нет
    with pytest.raises(api.RequestError, match="числом"):
        solve(2, MIN, r="x")

def test_supercritical_minimum_cost():
    _, body = solve(6, MIN, after="1")
    d = body["damping"]
    assert d["regime"] == "сверхкритический" and d["K"] == 3 and d["rho"] == pytest.approx(0.0)
    assert d["cost"] == pytest.approx(0.8 ** 2 / 6 / 3, rel=1e-9)          # J*/K
    assert d["exact"] and d["state_norm_at_T"] < 1e-25 and failed(body) == []
    assert set(d["count"]) == {3}
    assert all(f["t"][-1] <= 6.0 + 1e-9 for f in body["fronts"])           # после T разрывов нет

def test_subcritical_reports_unreachable_goal_with_minimal_residual():
    _, body = solve(1.5, MIN, after="1")
    d = body["damping"]
    expected = 0.8 ** 2 * (2 - 1.5) ** 3 / 24                               # κ²a³(2t0 − T)³/24
    assert d["regime"] == "докритический" and d["exact"] is False and d["cost_min"] is None
    assert d["state_norm_at_T"] == pytest.approx(expected, rel=1e-7)
    assert d["residual"] == pytest.approx(expected, rel=1e-9) and d["residual_min"] == pytest.approx(expected, rel=1e-9)
    assert d["state_norm_at_end"] == pytest.approx(expected, rel=1e-6)      # без управления норма сохраняется
    assert any("невозможно" in n for n in d["notes"])
    assert failed(body) == []
    names = [c["name"] for c in body["checks"]]
    assert "остаток Φ(T)" in names and "наименьший возможный остаток Φ_min(T)" in names

def test_subcritical_exact_case_is_recognised():
    preset = PRESETS["Успокоение: волна влево гасится за 0.6·t₀"]
    _, body = api.solve(dict(preset, nx=41, nt=41, fd_nodes=60))
    d = body["damping"]
    assert d["regime"] == "докритический" and d["exact"] is True
    assert d["minimal_time"] == pytest.approx(0.6, abs=2e-3)
    assert d["state_norm_at_T"] < 1e-20 and failed(body) == [] and body["sources"] == []

def test_continuous_designer():
    _, body = solve(3.3, CONT, eps="0.1", after="0.7")
    d = body["damping"]
    assert d["mode"] == "непрерывное" and d["exact"] is True
    assert body["sources"] == [] and body["norm"]["energy"] is not None     # скачков нет — энергия конечна
    assert d["cost_min"] < d["cost"] < 1.4 * d["cost_min"]
    assert body["boundary"]["u"][0] == pytest.approx(0.8)                   # u(0) = y0(0)
    assert failed(body) == []

def test_continuous_at_critical_time_is_the_consistent_control():
    _, body = solve(2, CONT, after="0.5")
    d = body["damping"]
    assert d["r"] == pytest.approx(0.8) and body["sources"] == []
    assert d["cost"] == pytest.approx(4 * d["cost_min"], rel=1e-9)
    assert any("единственно" in n for n in d["notes"])

def test_continuous_refuses_subcritical_time_and_warns_about_right_end():
    with pytest.raises(api.RequestError, match="однозначно"):
        solve(1.5, CONT)
    with pytest.raises(api.RequestError, match="положительной"):
        solve(3, CONT, eps="0")
    _, body = api.solve(dict(BASE, y0="0.8*(1 - x/l) + 0.3*x/l", T="3", source=CONT, options={"eps": "0.1"}))
    assert any("y₀(l) ≠ 0" in n for n in body["damping"]["notes"])
    assert body["damping"]["exact"] is True                                 # покой всё равно достигается

def test_two_passes_then_rest():
    _, body = solve(5, TWO, variant="оптимальное h", after="0")
    d = body["damping"]
    assert d["target_time"] == pytest.approx(2.0) and body["params"]["T"] == pytest.approx(5.0)
    assert d["cost"] == pytest.approx(0.8 ** 2 / 6, rel=1e-9) and d["state_norm_at_end"] < 1e-25
    _, body2 = solve(5, TWO, variant="согласованное (непрерывное)", after="0")
    assert body2["damping"]["cost"] == pytest.approx(4 * d["cost"], rel=1e-9) and body2["sources"] == []
    with pytest.raises(api.RequestError, match="2·t₀"):
        solve(1.5, TWO)

def test_extra_checks_are_registered_once_and_cover_formula_controls():
    names = [f.__name__ for f in diagnostics.EXTRA_CHECKS]
    assert names.count("linking_consistency") == 1 and names.count("damping_goal") == 1
    pr = Problem.from_text(l=1.0, a=1.0, T=2.7, y0="0.5*sin(pi*x/l)", y1="0.2", u="0.3*sin(4*t)")
    rows = [c for c in diagnostics.run_all(pr, pr.solve(), fd_nodes=0) if c.group == "Условие связи"]
    assert len(rows) == 1 and rows[0].ok is True and rows[0].value < 1e-9

def test_goal_check_fails_for_a_control_that_misses_the_goal():
    """Если управление помечено как точное, а покоя нет, проверка это показывает."""
    pr = Problem.from_text(l=1.0, a=1.0, T=2.0, y0="0.8*(1 - x/l)", y1="0", u="0.1*t")
    pr.control.damping = {"target_time": 2.0, "exact": True, "initial_norm": 0.2, "residual": 0.0, "residual_min": 0.0}
    rows = damping.checks.damping_goal(pr, pr.solve())
    assert rows[0].ok is False

@pytest.mark.parametrize("name", [n for n, p in PRESETS.items() if "source" in p])
def test_damping_presets_pass_all_checks(name):
    _, body = api.solve(dict(PRESETS[name], nx=41, nt=61, fd_nodes=80))
    assert failed(body) == []
    assert body["damping"] is not None
    json.dumps(body, allow_nan=False)

def test_formula_source_has_no_damping_report():
    _, body = api.solve(dict(BASE, T="2", source=api.FORMULA_SOURCE))
    assert body["damping"] is None

def test_regime_names():
    pr = Problem.from_text(l=1.0, a=1.0, T=1.0, y0="x", y1="0", u="0")
    for T, word in ((1.5, "докритический"), (2.0, "критический"), (2.0000001, "сверхкритический"), (7.0, "сверхкритический")):
        assert designers.regime(damping.BaseProfile(pr.params, pr.data, T))[0] == word
