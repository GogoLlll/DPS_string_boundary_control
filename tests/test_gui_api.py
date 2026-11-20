import csv
import io
import json

import numpy as np
import pytest

from wavesys import design
from wavesys.controls import FunctionControl
from wavesys.gui import api
from wavesys.gui.serialize import arr, num
from wavesys.presets import DEFAULT, PRESETS

BASE = dict(l="1", a="1", T="2", y0="0.5*sin(pi*x/l)", y1="0", u="0.1*t", nx=41, nt=31, fd_nodes=40)

def test_serialize_replaces_non_finite():
    assert num(float("nan")) is None and num(float("inf")) is None and num(2) == 2.0 and num(None) is None
    assert arr(np.array([1.0, np.nan, np.inf])) == [1.0, None, None]
    assert arr(np.array([[1.0, 2.0], [np.nan, 4.0]])) == [[1.0, 2.0], [None, 4.0]]
    json.dumps(arr(np.array([np.nan])), allow_nan=False)

def test_meta_lists_presets_and_formula_source():
    m = api.meta()
    assert m["default"] == DEFAULT and set(m["presets"]) == set(PRESETS)
    assert m["sources"][0]["name"] == api.FORMULA_SOURCE
    json.dumps(m, allow_nan=False)

def test_solve_returns_valid_json_with_expected_shapes():
    r, body = api.solve(BASE)
    text = json.dumps(body, allow_nan=False)                 # в ответе не должно быть NaN
    back = json.loads(text)
    assert len(back["grid"]["t"]) == 31 and len(back["grid"]["x"]) == 41
    assert len(back["grid"]["Y"]) == 31 and len(back["grid"]["Y"][0]) == 41
    assert back["params"] == {"l": 1.0, "a": 1.0, "t0": 1.0, "T": 2.0}
    assert back["stats"]["checks_passed"] == back["stats"]["checks_total"] > 0
    assert {"group", "name", "value", "limit", "ok", "note"} <= set(back["checks"][0])
    assert back["norm"]["energy"] is not None and back["sources"] == []

def test_solve_accepts_decimal_comma_and_auto_frames():
    payload = dict(BASE, l="1,5", nt="auto")
    r, body = api.solve(payload)
    assert body["params"]["l"] == 1.5
    assert len(body["grid"]["t"]) == api.auto_frames(2.0, 1.5)

def test_auto_frames_bounds():
    assert api.auto_frames(0.1, 1.0) == 121
    assert api.auto_frames(2.0, 1.0) == 121
    assert api.auto_frames(5.0, 1.0) == 301
    assert api.auto_frames(1000.0, 1.0) == 721

def test_discontinuous_problem_reports_fronts_and_no_energy():
    r, body = api.solve(dict(BASE, y0="0.5*(1 - x/l)", u="0"))
    assert body["norm"]["energy"] is None
    (src,) = body["sources"]
    assert src["kind"] == "left" and src["size"] == pytest.approx(-0.5)
    (front,) = body["fronts"]
    assert front["t"] == [0.0, 1.0, 2.0] and front["x"] == [0.0, 1.0, 0.0] and front["size"] == pytest.approx(-0.5)
    assert "скачки" in body["stats"]["verdict"]

@pytest.mark.parametrize("change, fragment", [
    (dict(l="abc"), "числом"), (dict(l=""), "Не задано"), (dict(a="-1"), "положительным"), (dict(T="0"), "положительным"),
    (dict(y0="sin("), "разобрать"), (dict(u="zzz(t)"), "неизвестная функция"), (dict(y1="t"), "неизвестное имя"),
    (dict(nx="5"), "меньше"), (dict(T="1e9"), "Слишком большое время"), (dict(source="нет такого"), "Неизвестный источник"),
    (dict(l="inf"), "конечным"),
])
def test_bad_input_gives_readable_error(change, fragment):
    with pytest.raises(api.RequestError, match=fragment):
        api.solve(dict(BASE, **change))

def test_point_matches_solution_and_explains_it():
    r, _ = api.solve(BASE)
    info = api.point(r, {"t": "1.3", "x": "0.4"})
    assert info["y"] == pytest.approx(float(r.solution.y(1.3, 0.4)), abs=1e-14)
    assert info["y"] == pytest.approx(0.5 * (info["A"] + info["B"]))
    assert info["kA"] == 1 and info["kB"] == 1 and info["region"] is None
    assert all(m["t"] <= 1.3 for m in info["control_A"] + info["control_B"])
    assert api.point(r, {"t": "0.2", "x": "0.5"})["region"] == "I"
    with pytest.raises(api.RequestError):
        api.point(r, {"t": "1", "x": "7"})                  # точка вне струны
    json.dumps(info, allow_nan=False)

def test_sections_are_exact_and_validated():
    r, _ = api.solve(BASE)
    s = api.sections(r, {"times": "0.5, 1.5", "points": "0.25"})
    x = np.array(s["x"])
    assert np.allclose(s["shapes"][1]["y"], r.solution.y(np.full(x.shape, 1.5), x))
    t = np.array(s["t"])
    assert np.allclose(s["traces"][0]["y"], r.solution.y(t, np.full(t.shape, 0.25)))
    for bad in ({"times": "", "points": "0.5"}, {"times": "9", "points": "0.5"}, {"times": "1", "points": "a"},
                {"times": "0,0.1,0.2,0.3,0.4,0.5", "points": "0.5"}):
        with pytest.raises(api.RequestError):
            api.sections(r, bad)

def test_table_and_csv_agree():
    r, _ = api.solve(BASE)
    t = api.table(r, {"rows": "4", "cols": "3"})
    rows = list(csv.reader(io.StringIO(api.table_csv(r, {"rows": "4", "cols": "3"}))))
    assert len(rows) == 5 and len(rows[0]) == 4
    assert np.allclose([[float(v) for v in row[1:]] for row in rows[1:]], t["Y"], atol=1e-9)

def test_registered_designer_is_offered_and_used():
    class Half(design.ControlDesigner):
        name = "тест: половина"
        description = "u = k * y0(0)"
        options = {"k": ("множитель", "0.5")}

        def build(self, params, data, horizon, k="0.5", **_):
            level = float(k) * float(np.asarray(data.y0(np.array([0.0]))).ravel()[0])
            return FunctionControl(lambda t: np.full(np.shape(t), level), label=f"{k}·y0(0)")

    design.register(Half())
    try:
        names = [s["name"] for s in api.meta()["sources"]]
        assert "тест: половина" in names
        r, body = api.solve(dict(BASE, y0="0.8*(1 - x/l)", source="тест: половина", options={"k": "0.25"}))
        assert body["boundary"]["u"][0] == pytest.approx(0.2)
        assert body["texts"]["u"] == "0.25·y0(0)"
    finally:
        design.REGISTRY.pop("тест: половина", None)
