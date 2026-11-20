import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from wavesys.characteristics import Source, front_position
from wavesys.params import StringParams

JS = Path(__file__).resolve().parent.parent / "wavesys" / "gui" / "static" / "js"
node = shutil.which("node")
if node is None and os.environ.get("REQUIRE_NODE"):      # в CI отсутствие node — ошибка, а не пропуск
    raise RuntimeError("REQUIRE_NODE задан, но node не найден")
pytestmark = pytest.mark.skipif(node is None, reason="node не установлен")


def run_js(code):
    script = f'import * as F from "{(JS / "fronts.js").as_posix()}";\nimport * as N from "{(JS / "format.js").as_posix()}";\n{code}'
    out = subprocess.run([node, "--input-type=module", "-e", script], capture_output=True, text=True, timeout=30)
    assert out.returncode == 0, out.stderr
    return json.loads(out.stdout)

def test_all_scripts_parse():
    for path in sorted(JS.rglob("*.js")):
        out = subprocess.run([node, "--check", "--input-type=module"], input=path.read_text(encoding="utf-8"),
                             capture_output=True, text=True, timeout=30)
        assert out.returncode == 0, f"{path.name}: {out.stderr}"

def test_break_between_points_extends_both_sides_to_the_cut():
    r = run_js("console.log(JSON.stringify(F.breakAt([0,1,2,3,4],[10,11,12,13,14],[1.5])))")
    assert r == {"x": [0, 1, 1.5, 1.5, 1.5, 2, 3, 4], "y": [10, 11, 11, None, 12, 12, 13, 14]}

def test_point_on_the_cut_is_replaced_by_one_sided_values():
    r = run_js("console.log(JSON.stringify(F.breakAt([0,1,2,3,4],[10,11,12,13,14],[2])))")
    assert r == {"x": [0, 1, 2, 2, 2, 3, 4], "y": [10, 11, 11, None, 13, 13, 14]}

def test_cuts_at_ends_and_outside():
    r = run_js("console.log(JSON.stringify([F.breakAt([0,1,2],[5,6,7],[0,2]), F.breakAt([0,1,2],[5,6,7],[null,9,-3])]))")
    assert r[0]["y"][0] is None and r[0]["y"][-1] is None and 6 in r[0]["y"]
    assert r[1] == {"x": [0, 1, 2], "y": [5, 6, 7]}

def test_front_position_agrees_with_python():
    p = StringParams(2.0, 4.0)
    times = [0.0, 0.1, 0.37, 0.5, 0.81, 1.0, 1.9, 3.3]
    for kind, start in (("left", 0.0), ("right", 0.0), ("left", 0.3)):
        js = run_js(f'const p={{l:2,a:4,t0:0.5}};console.log(JSON.stringify({json.dumps(times)}.map(t=>F.frontPosition(p,{{kind:"{kind}",time:{start}}},t))))')
        py = front_position(p, Source(kind, start, 1.0), np.array(times))
        for a, b in zip(js, py):
            assert (a is None and np.isnan(b)) or a == pytest.approx(float(b), abs=1e-12)

def test_crossing_times_are_where_the_front_passes_the_point():
    p = StringParams(2.0, 4.0)
    for kind in ("left", "right"):
        js = run_js(f'console.log(JSON.stringify(F.crossingTimes({{l:2,a:4,t0:0.5}},{{kind:"{kind}",time:0.1}},0.7,3)))')
        assert len(js) >= 5
        pos = front_position(p, Source(kind, 0.1, 1.0), np.array(js))
        assert np.allclose(pos, 0.7, atol=1e-12)

def test_profile_jumps_follow_the_reflection_law():
    js = run_js('console.log(JSON.stringify([F.profileJumps({l:2,a:4,t0:0.5},{kind:"left",time:0.2},2.2), F.profileJumps({l:2,a:4,t0:0.5},{kind:"right",time:0},1.5)]))')
    assert np.allclose(js[0]["B"], [0.7, 1.7]) and np.allclose(js[0]["A"], [1.2, 2.2])     # B в τ + t0, A в τ + 2 t0, далее через 2 t0
    assert np.allclose(js[1]["A"], [0.5, 1.5]) and np.allclose(js[1]["B"], [1.0])

def test_number_formatting():
    r = run_js('console.log(JSON.stringify([N.fmt(0), N.fmt(0.123456), N.fmt(-0.5), N.fmt(1234.5678), N.fmt(1.5e-7), N.fmt(null), N.sci(-2.5e-16), N.inT0(1.5, 0.75)]))')
    assert r == ["0", "0.1235", "−0.5", "1235", "1.50·10⁻⁷", "—", "−2.50·10⁻¹⁶", "2 t₀"]
