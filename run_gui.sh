#!/bin/bash
# Запуск интерфейса. Дополнительные параметры: --port N, --no-browser
cd "$(dirname "$0")" && exec .venv/bin/python -m wavesys.gui "$@"
