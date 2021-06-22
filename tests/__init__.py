# -*- coding: utf-8 -*-
""" Tests for the plugin.

Includes both tests written in unittest style (test_cli.py) and tests written
in pytest style (test_calculations.py).
"""
from pathlib import Path
import os

TEST_DIR = Path(__file__).resolve().parent
DATA_DIR = TEST_DIR / 'data'
