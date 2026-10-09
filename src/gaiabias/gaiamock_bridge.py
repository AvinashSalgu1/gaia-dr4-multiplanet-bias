"""Locate and import gaiamock, which is not pip-installable.

Search order for the gaiamock directory:
  1. the GAIAMOCK_DIR environment variable;
  2. ../gaiamock next to this repository (the recommended layout:
     <parent>/gaia-dr4-multiplanet-bias and <parent>/gaiamock).
"""
from __future__ import annotations

import importlib
import os
import sys
from functools import lru_cache
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def gaiamock_dir() -> Path:
    env = os.environ.get("GAIAMOCK_DIR")
    path = Path(env).expanduser() if env else REPO_ROOT.parent / "gaiamock"
    return path.resolve()


def available() -> bool:
    """True if gaiamock imports and its compiled Kepler solver exists."""
    try:
        gaiamock()
        c_funcs()
        return True
    except Exception:
        return False


@lru_cache(maxsize=1)
def gaiamock():
    d = gaiamock_dir()
    if not (d / "gaiamock.py").exists():
        raise ImportError(f"gaiamock.py not found in {d}. Set GAIAMOCK_DIR.")
    if str(d) not in sys.path:
        sys.path.insert(0, str(d))
    return importlib.import_module("gaiamock")


@lru_cache(maxsize=1)
def c_funcs():
    return gaiamock().read_in_C_functions()


def n_scan_files() -> int:
    d = gaiamock_dir() / "healpix_scans"
    return len(list(d.glob("healpix_64_*.fits"))) if d.exists() else 0
