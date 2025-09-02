#!/usr/bin/env python
# coding: utf-8

import shutil

def terminal_width(default: int = 80) -> int:
    """Return current terminal width, fallback to default."""
    try:
        return shutil.get_terminal_size().columns
    except Exception:
        return default

def print_full_line(char: str = "=", width: int | None = None) -> None:
    """Print a line composed of `char` filling the terminal width."""
    w = width or terminal_width()
    print(char * max(0, w))
