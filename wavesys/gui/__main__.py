import argparse
import subprocess
import sys
import webbrowser
from pathlib import Path

from .server import make_server

def open_browser(url):
    try:
        if "microsoft" in Path("/proc/version").read_text().lower():
            subprocess.Popen(["explorer.exe", url], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return
    except OSError:
        pass
    try:
        webbrowser.open(url)
    except Exception:
        pass

def main(argv=None):
    parser = argparse.ArgumentParser(description="Интерфейс wavesys: волновое уравнение с граничным управлением")
    parser.add_argument("--port", type=int, default=8765, help="порт (по умолчанию 8765; если занят — следующий свободный)")
    parser.add_argument("--no-browser", action="store_true", help="не открывать браузер")
    args = parser.parse_args(argv)
    server = make_server(args.port)
    url = f"http://localhost:{server.server_address[1]}/"
    if not args.no_browser:
        open_browser(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0

if __name__ == "__main__":
    sys.exit(main())
