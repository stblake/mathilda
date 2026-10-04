#!/usr/bin/env python3
"""Documentation coverage report for the builtin-doc overhaul.

A page is "EllipticF-grade" when all three hand-authored layers exist:

  * a worked-example + Notes overlay  (site/overlays/<Name>.md)
  * a source-grounded implementation note  (site/impl/<Name>.md)
  * at least one verified example on the generated page
    (the page does not say "No verified examples yet")

This reads the *generated* tree under site/docs/documentation/ (the state of
the last `make docs`) plus the overlay/impl source dirs, and prints a
per-category table of (pages, overlay, impl, examples, EllipticF-grade). It is
the progress metric for BUILTIN_DOCUMENTATION_OVERHAUL.md.

Usage:
  python3 site/coverage_report.py                 # full per-category table
  python3 site/coverage_report.py <category-slug> # one category, listing gaps
  python3 site/coverage_report.py --missing       # every function missing a layer
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DOC = ROOT / "docs" / "documentation"
OVERLAYS = ROOT / "overlays"
IMPL = ROOT / "impl"

NO_EXAMPLES = "No verified examples yet"


def scan():
    """Return {category: [(name, has_overlay, has_impl, has_examples), ...]}."""
    cats = {}
    for page in sorted(DOC.glob("*/*.md")):
        if page.name == "index.md":
            continue
        name = page.stem
        cat = page.parent.name
        text = page.read_text(errors="replace")
        row = (
            name,
            (OVERLAYS / f"{name}.md").exists(),
            (IMPL / f"{name}.md").exists(),
            NO_EXAMPLES not in text,
        )
        cats.setdefault(cat, []).append(row)
    return cats


def graded(row):
    _, ov, im, ex = row
    return ov and im and ex


def print_table(cats):
    hdr = f"{'category':<26}{'pages':>7}{'overlay':>9}{'impl':>7}{'examples':>10}{'graded':>8}{'%':>6}"
    print(hdr)
    print("-" * len(hdr))
    tot = [0, 0, 0, 0, 0]
    for cat in sorted(cats, key=lambda c: -len(cats[c])):
        rows = cats[cat]
        n = len(rows)
        ov = sum(r[1] for r in rows)
        im = sum(r[2] for r in rows)
        ex = sum(r[3] for r in rows)
        gr = sum(graded(r) for r in rows)
        tot = [tot[0] + n, tot[1] + ov, tot[2] + im, tot[3] + ex, tot[4] + gr]
        pct = 100 * gr // n if n else 0
        print(f"{cat:<26}{n:>7}{ov:>9}{im:>7}{ex:>10}{gr:>8}{pct:>5}%")
    print("-" * len(hdr))
    n = tot[0]
    pct = 100 * tot[4] // n if n else 0
    print(f"{'TOTAL':<26}{n:>7}{tot[1]:>9}{tot[2]:>7}{tot[3]:>10}{tot[4]:>8}{pct:>5}%")


def print_category(cats, slug):
    rows = cats.get(slug)
    if not rows:
        print(f"no such category: {slug}", file=sys.stderr)
        print("categories:", ", ".join(sorted(cats)), file=sys.stderr)
        sys.exit(2)
    print(f"# {slug} ({len(rows)} functions)\n")
    print(f"{'function':<30}{'overlay':>9}{'impl':>7}{'examples':>10}{'graded':>8}")
    for r in sorted(rows):
        name, ov, im, ex = r
        mark = lambda b: " yes" if b else "  --"
        g = " yes" if graded(r) else "  --"
        print(f"{name:<30}{mark(ov):>9}{mark(im):>7}{mark(ex):>10}{g:>8}")
    gaps = [r[0] for r in sorted(rows) if not graded(r)]
    print(f"\n{len(gaps)} not yet EllipticF-grade:")
    print("  " + " ".join(gaps) if gaps else "  (all graded)")


def print_missing(cats):
    for cat in sorted(cats):
        for r in sorted(cats[cat]):
            if graded(r):
                continue
            name, ov, im, ex = r
            need = []
            if not ov:
                need.append("overlay")
            if not im:
                need.append("impl")
            if not ex:
                need.append("examples")
            print(f"{cat}/{name}: needs {', '.join(need)}")


def main():
    cats = scan()
    args = sys.argv[1:]
    if not args:
        print_table(cats)
    elif args[0] == "--missing":
        print_missing(cats)
    else:
        print_category(cats, args[0])


if __name__ == "__main__":
    main()
