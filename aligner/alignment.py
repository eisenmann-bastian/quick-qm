#!/usr/bin/env python3
"""Compatibility entry point for deriving axes from a selector CSV."""

import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from aligner.axes import main

if __name__ == "__main__":
    sys.exit(main())
