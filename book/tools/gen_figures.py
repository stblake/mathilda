#!/usr/bin/env python3
"""Generate the book's plot cells: a REPL session that ends in a graphic.

Each ``book/figures/<path>.m`` is one Mathilda session whose final line binds a
plot to ``fig`` (``fig = Plot[...]``). This tool runs the session through the
real ``./Mathilda`` and writes two files under ``book/generated/figures/``:

  <path>.in    the In[k]/Out[k] transcript, EXACTLY as build_examples renders
               it, except the plot line is shown as its bare command (the
               ``fig = `` is stripped, so the reader sees ``Plot[...]``) and its
               output is the label ``Out[k]= `` with no text -- the graphic goes
               there instead.
  <path>.pdf   that plot, Export-ed to a vector PDF by the binary itself.

The book includes both with ``\\plotcell{<path>}`` (book/mathilda.sty): the
transcript in a REPL box, then the PDF centred right beneath ``Out[k]=``. So the
figure is the genuine output of a shown command, produced the same way every
other transcript in the book is -- nothing hand-drawn, nothing hand-typed.

PDF export is the dependency-free vector path (no window/display), so this runs
headless in CI alongside the rest of the build.
"""
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]          # repo root
sys.path.insert(0, str(ROOT / "site"))
from verify_tutorial import run_session, MATHILDA    # noqa: E402 (reuse the driver)
sys.path.insert(0, str(ROOT / "book" / "tools"))
from build_examples import read_inputs               # noqa: E402 (same input parser)

FIGSRC = ROOT / "book" / "figures"
OUT = ROOT / "book" / "generated" / "figures"

FIG_ASSIGN = re.compile(r"^\s*fig\s*=\s*")
# Optional first-line directive: choose the output format. Vector PDF is the
# default (crisp, tiny for line plots); PNG is for the dense colour plots
# (DensityPlot/ComplexPlot/ArrayPlot -- megabytes as vectors) and all 3D
# (Graphics3D has no vector form). e.g.  # format: png
FORMAT_RE = re.compile(r"^#\s*format:\s*(pdf|png)\s*$", re.I)


def figure_format(path):
    for ln in path.read_text().splitlines():
        m = FORMAT_RE.match(ln.strip())
        if m:
            return m.group(1).lower()
    return "pdf"


def render_cell(inputs, outs, plot_idx):
    """One In[k]/Out[k] block per input line. The plot line shows its bare
    command and an empty ``Out[k]= `` label (the image is placed there by the
    macro); other lines follow build_examples' convention (Out omitted when the
    result is empty or Null)."""
    blocks = []
    for k, inp in enumerate(inputs, 1):
        i = k - 1
        if i == plot_idx:
            cmd = FIG_ASSIGN.sub("", inp)
            blocks.append(f"In[{k}]:= {cmd}\nOut[{k}]= ")
        else:
            out = (outs[i].strip() if i < len(outs) else "")
            line = f"In[{k}]:= {inp}"
            if out not in ("", "Null"):
                line += f"\nOut[{k}]= {out}"
            blocks.append(line)
    return "\n\n".join(blocks) + "\n"


def main():
    scripts = sorted(FIGSRC.rglob("*.m")) if FIGSRC.exists() else []
    if not scripts:
        print("gen_figures: no book/figures/**/*.m found")
        return
    if not MATHILDA.exists():
        sys.exit(f"error: {MATHILDA} not found -- run `make` in the repo root first.")

    # Belt and braces: PDF export opens no window, but make sure a stray
    # auto-display can never try to either.
    os.environ["MATHILDA_NO_GRAPHICS_WINDOW"] = "1"

    n_ok = 0
    for m in scripts:
        rel = m.relative_to(FIGSRC).with_suffix("")     # e.g. 02-introduction/ndsolve-solution
        ext = figure_format(m)                           # "pdf" (default) or "png"
        img = OUT / f"{rel}.{ext}"
        dst = OUT / f"{rel}.in"
        img.parent.mkdir(parents=True, exist_ok=True)
        # Drop any stale sibling in the other format so \plotcell never picks it.
        other = OUT / f"{rel}.{'png' if ext == 'pdf' else 'pdf'}"
        if other.exists():
            other.unlink()

        inputs = read_inputs(m)
        plot_idx = next((i for i, ln in enumerate(inputs) if FIG_ASSIGN.match(ln)), None)
        if plot_idx is None:
            sys.exit(f"gen_figures: {m.relative_to(ROOT)} has no `fig = <plot>` line")

        # Run the session and Export `fig`; capture In/Out for the shown lines.
        outs = run_session(inputs + [f'Export["{img}", fig]'])
        dst.write_text(render_cell(inputs, outs, plot_idx))

        if img.exists() and img.stat().st_size > 0:
            print(f"gen_figures: {m.relative_to(ROOT)} -> {img.relative_to(ROOT)} "
                  f"({img.stat().st_size} bytes) + {dst.name}")
            n_ok += 1
        else:
            sys.exit(f"gen_figures: FAILED to produce {img} from {m.relative_to(ROOT)} "
                     f"(is `fig` bound to a Graphics object? PNG needs a GUI session.)")
    print(f"gen_figures: wrote {n_ok} plot cell(s) -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
