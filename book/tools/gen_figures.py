#!/usr/bin/env python3
"""Generate the book's figures by running figure scripts through ./Mathilda and
Export-ing the resulting plot to a PDF under ``book/generated/figures/``.

Design contract (mirrors build_examples.py): a figure is never a hand-drawn
asset. Each ``book/figures/<path>.m`` is ONE Mathilda session that builds a plot
and binds it to the variable ``fig``; this tool appends
``Export["<abs>/generated/figures/<path>.pdf", fig]`` and runs it, so the figure
is produced by the real binary exactly like every ``In[]/Out[]`` transcript. PDF
export is the dependency-free vector path (no window, no display), so this runs
headless in CI alongside the rest of the build.

The book includes each result with ``\\bookfigure{<path>}{caption}`` (see
book/mathilda.sty), which falls back to a visible placeholder when the figure has
not been generated yet, so the document always compiles.
"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]          # repo root
sys.path.insert(0, str(ROOT / "site"))
from verify_tutorial import run_session, MATHILDA    # noqa: E402 (reuse the driver)

FIGSRC = ROOT / "book" / "figures"
OUT = ROOT / "book" / "generated" / "figures"


def figure_scripts():
    return sorted(FIGSRC.rglob("*.m")) if FIGSRC.exists() else []


def source_lines(path):
    """Non-blank, non-comment input lines, one expression per line."""
    return [ln for ln in path.read_text().splitlines()
            if ln.strip() and not ln.lstrip().startswith("#")]


def main():
    scripts = figure_scripts()
    if not scripts:
        print("gen_figures: no book/figures/**/*.m found")
        return
    if not MATHILDA.exists():
        sys.exit(f"error: {MATHILDA} not found -- run `make` in the repo root first.")

    # Belt and braces: even though PDF export opens no window, make sure a stray
    # auto-display can never try to.
    os.environ["MATHILDA_NO_GRAPHICS_WINDOW"] = "1"

    n_ok = 0
    for m in scripts:
        rel = m.relative_to(FIGSRC).with_suffix("")     # e.g. 02-introduction/ndsolve-solution
        pdf = OUT / f"{rel}.pdf"
        pdf.parent.mkdir(parents=True, exist_ok=True)
        lines = source_lines(m)
        lines.append(f'Export["{pdf}", fig]')
        run_session(lines)
        if pdf.exists() and pdf.stat().st_size > 0:
            print(f"gen_figures: {m.relative_to(ROOT)} -> {pdf.relative_to(ROOT)} "
                  f"({pdf.stat().st_size} bytes)")
            n_ok += 1
        else:
            sys.exit(f"gen_figures: FAILED to produce {pdf} from {m.relative_to(ROOT)} "
                     f"(is `fig` bound to a Graphics object?)")
    print(f"gen_figures: wrote {n_ok} figure(s) -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
