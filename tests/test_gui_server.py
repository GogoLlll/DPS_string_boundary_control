import json
import threading
import urllib.error
import urllib.request

import pytest

from wavesys.gui.server import make_server

@pytest.fixture(scope="module")
def base_url():
    server = make_server(port=0)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{server.server_address[1]}"
    server.shutdown()
    server.server_close()

def get(url):
    with urllib.request.urlopen(url, timeout=20) as resp:
        return resp.status, resp.headers, resp.read()

def post(url, payload):
    req = urllib.request.Request(url, data=json.dumps(payload).encode(), headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return resp.status, json.loads(resp.read())

def error_of(call):
    with pytest.raises(urllib.error.HTTPError) as exc:
        call()
    return exc.value.code, json.loads(exc.value.read())["error"]

def test_page_and_static_files(base_url):
    status, headers, body = get(base_url + "/")
    assert status == 200 and "text/html" in headers["Content-Type"] and "problem-form".encode() in body
    for path, kind in (("/style.css", "text/css"), ("/js/main.js", "javascript"), ("/js/views/overview.js", "javascript")):
        status, headers, body = get(base_url + path)
        assert status == 200 and kind in headers["Content-Type"] and len(body) > 100

def test_every_module_imported_by_the_page_exists(base_url):
    import re
    _, _, body = get(base_url + "/js/main.js")
    paths = re.findall(r'from "\./([\w/]+\.js)"', body.decode())
    assert len(paths) >= 10
    for rel in paths:
        assert get(f"{base_url}/js/{rel}")[0] == 200

def test_plotting_library_is_served_locally(base_url):
    status, _, body = get(base_url + "/vendor/plotly.min.js")
    assert status == 200 and len(body) > 1_000_000

def test_cannot_escape_static_directory(base_url):
    for path in ("/../server.py", "/js/../../server.py", "/%2e%2e/server.py", "/nothing.txt"):
        code, _ = error_of(lambda: get(base_url + path))
        assert code == 404

def test_requests_before_first_solve_are_refused(base_url):
    code, message = error_of(lambda: get(base_url + "/api/point?t=1&x=0.5"))
    assert code == 400 and "Сначала" in message

def test_solve_then_follow_up_requests(base_url):
    status, body = post(base_url + "/api/solve", dict(l=1, a=1, T=2, y0="0.5*sin(pi*x/l)", y1="0", u="0", nx=31, nt=21, fd_nodes=30))
    assert status == 200 and body["stats"]["checks_passed"] == body["stats"]["checks_total"]
    _, _, raw = get(base_url + "/api/point?t=0.5&x=0.5")
    assert json.loads(raw)["y"] == pytest.approx(0.0, abs=1e-12)          # 0.5 sin(π/2) cos(π/2)
    _, _, raw = get(base_url + "/api/table?rows=3&cols=3")
    assert json.loads(raw)["Y"][0][1] == pytest.approx(0.5)
    status, headers, raw = get(base_url + "/api/table.csv?rows=3&cols=3")
    assert "text/csv" in headers["Content-Type"] and "attachment" in headers["Content-Disposition"]
    assert raw.decode("utf-8-sig").splitlines()[0].startswith("t \\ x,0,")
    _, _, raw = get(base_url + "/api/sections?times=0,1&points=0.5")
    assert len(json.loads(raw)["shapes"]) == 2

def test_errors_are_reported_as_json(base_url):
    code, message = error_of(lambda: post(base_url + "/api/solve", dict(l=1, a=1, T=2, y0="oops(x)", y1="0", u="0")))
    assert code == 400 and "oops" in message
    code, message = error_of(lambda: post(base_url + "/api/nothing", {}))
    assert code == 404
    req = urllib.request.Request(base_url + "/api/solve", data=b"not json", headers={"Content-Type": "application/json"})
    code, message = error_of(lambda: urllib.request.urlopen(req, timeout=20))
    assert code == 400 and "JSON" in message
    req = urllib.request.Request(base_url + "/api/solve", data=b"[1, 2]", headers={"Content-Type": "application/json"})
    assert error_of(lambda: urllib.request.urlopen(req, timeout=20))[0] == 400

def test_port_fallback_when_busy():
    first = make_server(port=0)
    try:
        busy = first.server_address[1]
        second = make_server(port=busy)
        try:
            assert second.server_address[1] != busy
        finally:
            second.server_close()
    finally:
        first.server_close()
