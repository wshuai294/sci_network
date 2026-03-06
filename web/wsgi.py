#!/usr/bin/env python3
"""WSGI entry for production (e.g. gunicorn). Run from project root: gunicorn -w 4 -b 0.0.0.0:5050 "web.wsgi:app" """
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "web"))

from app import app
