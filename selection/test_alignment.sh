#!/bin/sh
set -eu
cd "$(dirname "$0")/.."
export PYTHONDONTWRITEBYTECODE=1
exec "${PYTHON:-python3}" -m unittest selection.test_alignment "$@"