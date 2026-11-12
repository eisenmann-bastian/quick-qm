#!/bin/sh
set -eu
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
exec "${PYTHON:-python3}" -m unittest aligner.test_selector aligner.test_alignment aligner.test_axis_config aligner.test_aligner aligner.test_torsion_optimization aligner.test_cross_repository "$@"
