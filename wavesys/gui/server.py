import json
import mimetypes
import threading
import traceback
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from . import api

STATIC = Path(__file__).resolve().parent / "static"

def _plotly_path():
    import plotly
    return Path(plotly.__file__).resolve().parent / "package_data" / "plotly.min.js"

class AppState:
    def __init__(self):
        self.lock = threading.Lock()
        self.results = None


class Handler(BaseHTTPRequestHandler):
    state: AppState = None
    server_version = "wavesys"

    def log_message(self, fmt, *args):
        pass

    def _send(self, code, body, content_type, extra=None):
        if isinstance(body, str):
            body = body.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(body)

    def _json(self, obj, code=200):
        self._send(code, json.dumps(obj, ensure_ascii=False, allow_nan=False), "application/json; charset=utf-8")

    def _error(self, code, message):
        self._json({"error": message}, code)

    def _results(self):
        with self.state.lock:
            r = self.state.results
        if r is None:
            raise api.RequestError("Сначала выполните расчёт")
        return r

    def do_GET(self):
        url = urlparse(self.path)
        query = {k: v[0] for k, v in parse_qs(url.query).items()}
        try:
            if url.path == "/api/meta":
                return self._json(api.meta())
            if url.path == "/api/point":
                return self._json(api.point(self._results(), query))
            if url.path == "/api/sections":
                return self._json(api.sections(self._results(), query))
            if url.path == "/api/table":
                return self._json(api.table(self._results(), query))
            if url.path == "/api/table.csv":
                text = api.table_csv(self._results(), query)
                return self._send(200, "﻿" + text, "text/csv; charset=utf-8",
                                  {"Content-Disposition": 'attachment; filename="values.csv"'})
            if url.path == "/vendor/plotly.min.js":
                return self._send(200, _plotly_path().read_bytes(), "application/javascript; charset=utf-8")
            return self._static(url.path)
        except api.RequestError as exc:
            self._error(400, str(exc))
        except Exception as exc:
            traceback.print_exc()
            self._error(500, f"Внутренняя ошибка: {exc}")

    def do_POST(self):
        url = urlparse(self.path)
        try:
            length = int(self.headers.get("Content-Length") or 0)
            if length > 1_000_000:
                return self._error(413, "Слишком большой запрос")
            try:
                payload = json.loads(self.rfile.read(length) or b"{}")
            except json.JSONDecodeError:
                return self._error(400, "Запрос не является JSON")
            if not isinstance(payload, dict):
                return self._error(400, "Ожидается объект JSON")
            if url.path == "/api/solve":
                r, body = api.solve(payload)
                with self.state.lock:
                    self.state.results = r
                return self._json(body)
            self._error(404, "Нет такого адреса")
        except api.RequestError as exc:
            self._error(400, str(exc))
        except Exception as exc:
            traceback.print_exc()
            self._error(500, f"Внутренняя ошибка: {exc}")

    def _static(self, path):
        rel = "index.html" if path in ("", "/") else path.lstrip("/")
        target = (STATIC / rel).resolve()
        if STATIC not in target.parents or not target.is_file():
            return self._error(404, "Нет такого адреса")
        ctype = mimetypes.guess_type(target.name)[0] or "application/octet-stream"
        if target.suffix in (".html", ".css", ".js"):
            ctype += "; charset=utf-8"
        self._send(200, target.read_bytes(), ctype)


def make_server(port=8765, host="127.0.0.1", tries=20):
    state = AppState()
    handler = type("BoundHandler", (Handler,), {"state": state})
    last = None
    for candidate in ([0] if port == 0 else range(port, port + tries)):
        try:
            return ThreadingHTTPServer((host, candidate), handler)
        except OSError as exc:
            last = exc
    raise OSError(f"Не удалось занять порт {port}–{port + tries - 1}: {last}")
