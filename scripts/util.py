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
    
def print_box(text: str, padding: int = 5, min_width: int = 60) -> None:
    """
    Print `text` centered in a nice double-line box using Unicode box-drawing chars.

    Parameters
    - text: the string to display (can contain newlines; each line will be boxed)
    - padding: spaces left/right of text inside the box
    - min_width: minimum inner width (excluding vertical borders)
    """
    # Box characters (double-line)
    TL = "╔"
    TR = "╗"
    BL = "╚"
    BR = "╝"
    H  = "═"
    V  = "║"
    # For inner separators (if multiple lines) you could use double-line horizontal with middle joints
    # but here we keep a continuous top/bottom and vertical sides.

    # Split input into lines and compute inner width
    lines = text.splitlines() or [""]
    max_line_len = max(len(line) for line in lines)
    inner_width = max(max_line_len + 2 * padding, min_width)

    # Top border
    print(f"{TL}{H * inner_width}{TR}")

    # Content lines (center each input line)
    for line in lines:
        # compute centered placement
        space_total = inner_width - len(line)
        left_space = space_total // 2
        right_space = space_total - left_space
        print(f"{V}{' ' * left_space}{line}{' ' * right_space}{V}")

    # Bottom border
    print(f"{BL}{H * inner_width}{BR}")
