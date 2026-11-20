import csv

import numpy as np
import pytest

from wavesys import design, results
from wavesys.controls import FunctionControl
from wavesys.problem import Problem

def test_compute_shapes(smooth):
    r = results.compute(smooth, nx=31, nt=21, fd_nodes=40, norm_points=9)
    assert r.Y.shape == (21, 31) and r.t.shape == (21,) and r.x.shape == (31,)
    assert r.u.shape == r.left.shape == r.right.shape == (21,)
    assert r.A.shape == r.B.shape == r.z.shape
    assert r.norm.shape == r.norm_t.shape == (9,)
    assert r.z[-1] == pytest.approx(smooth.horizon + smooth.params.t0)

def test_grid_boundaries_match_control(smooth):
    r = results.compute(smooth, nx=31, nt=21, fd_nodes=0, norm_points=5)
    assert np.allclose(r.Y[:, 0], r.u, atol=1e-13)
    assert np.allclose(r.Y[:, -1], 0.0, atol=1e-13)
    assert np.allclose(r.left, r.u, atol=1e-13)

def test_energy_only_for_continuous_solutions(smooth):
    assert results.compute(smooth, nx=11, nt=11, fd_nodes=0, norm_points=5).energy is not None
    jumpy = Problem.from_text(l=1.0, a=1.0, T=1.0, y0="0.5*(1 - x/l)", y1="0", u="0")
    assert results.compute(jumpy, nx=11, nt=11, fd_nodes=0, norm_points=5).energy is None

def test_value_table_and_csv(tmp_path, smooth):
    r = results.compute(smooth, nx=21, nt=21, fd_nodes=0, norm_points=5)
    t, x, Y = results.value_table(r, rows=5, cols=4)
    assert Y.shape == (5, 4)
    assert np.allclose(Y, r.solution.grid(t, x))
    path = tmp_path / "table.csv"
    results.export_csv(path, r, rows=5, cols=4)
    rows = list(csv.reader(open(path, encoding="utf-8")))
    assert len(rows) == 6 and len(rows[0]) == 5
    assert np.allclose([float(v) for v in rows[0][1:]], x)
    assert np.allclose([[float(v) for v in row[1:]] for row in rows[1:]], Y, atol=1e-9)
    assert np.allclose([float(row[0]) for row in rows[1:]], t)

def test_designer_registry_roundtrip(smooth):
    class Hold(design.ControlDesigner):
        name = "тест: удержание"
        description = "держит левый конец на начальном смещении"

        def build(self, params, data, horizon, **options):
            level = float(np.asarray(data.y0(np.array([0.0]))).ravel()[0])
            return FunctionControl(lambda t: np.full(np.shape(t), level), label="удержание")

    try:
        d = design.register(Hold())
        assert design.REGISTRY["тест: удержание"] is d
        control = d.build(smooth.params, smooth.data, smooth.horizon)
        pr = Problem(smooth.params, smooth.data, control, smooth.horizon)
        assert pr.solve().y(0.0, 0.0) == pytest.approx(0.3)
        assert np.allclose(pr.solve().left_boundary(np.linspace(0, 2, 9)), 0.3)
    finally:
        design.REGISTRY.pop("тест: удержание", None)
