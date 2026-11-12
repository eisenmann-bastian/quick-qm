#!/usr/bin/env python3
"""Compatibility entry point for deriving axes from a selector CSV."""

import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from alignment.axes import (
    AlignmentCalculator, AlignmentCsv, AlignmentSdf, GROUP_COLUMNS, INPUT_COLUMNS,
    OUTPUT_COLUMNS, SDF_AXES_PROPERTY, SDF_SOURCE_PROPERTY, main,
)

if __name__ == "__main__":
    sys.exit(main())
