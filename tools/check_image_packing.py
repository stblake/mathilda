#!/usr/bin/env python3
"""check_image_packing.py -- does every image head hand back a PACKED buffer?

WHY THIS EXISTS. Three separate times in this subsystem an image operation was several times
slower than its equivalent elsewhere, and each time the algorithm was fine and the marshalling was
the cost:

  * `image_load` walked an NDArray element by element through `ndt_get`, found when ImagePad
    measured 0.57 ms against NumPy's 0.050.
  * `image3d_load` was still walking after `image_load` was fixed, found when the volumetric pad
    measured 6.5x NumPy while the planar one had just reached 2.2x.
  * `bit_image_from_mask` built 262144 `Expr` integers in nested `List`s, found when
    LocalAdaptiveBinarize measured 4.2x scikit-image -- and fixing it made global `Binarize`, an
    unrelated head, 23x faster.

In all three the ANSWERS WERE CORRECT and only the representation was wrong, so no test in the
suite could have caught them; a benchmark caught each one, by accident, one at a time. That is the
argument for a gate: the property "an image-returning head hands back a packed buffer" is
mechanically checkable, and checking it is much cheaper than benchmarking every head and noticing.

WHAT IT DOES. Heads are read out of the image sources rather than listed here, so a new one is
covered the day it is registered. For each head it DISCOVERS by trial which call shapes are
accepted -- the same approach `nd_fastpath_sweep.py` takes, and for the same reason: a curated list
of call shapes only covers the heads someone remembered. Any shape that returns an image is then
asked one question: is the stored data an `NDArray`, or is it nested `List`s?

It RATCHETS. `KNOWN_UNPACKED` carries the heads that are known to return nested data, so the gate
fails on a head that has newly fallen out and reports one that has been fixed. A gate that fails on
its whole standing backlog from the day it lands stops being read within a week.
"""

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "Mathilda"
SOURCES = ["src/image.c", "src/imagefilter.c", "src/imagegeom.c"]

# Heads known to return nested data. Each line is a standing debt, not a permission: delete the line
# when the head is fixed, and the gate will tell you to.
KNOWN_UNPACKED = {}

# Heads that legitimately do not return an image, so there is nothing to pack. Listed rather than
# inferred so that a head which STOPS returning an image is noticed.
NOT_IMAGE_RETURNING = {
    "ImageData",
    "ImageDimensions",
    "ImageChannels",
    "ImageType",
    "ImageQ",
    "Image3DQ",
    "ImageCorners",
    "FindThreshold",
    "ImageLevels",
    "GaussianMatrix",
    "BoxMatrix",
    "MorphologicalComponents",  # a matrix of integer labels, not an image
}

# Call shapes to try, as Mathilda source with `IMG` standing for the test image. Ordered cheapest
# first; the first shape that yields an image is the one reported.
SHAPES = [
    "{H}[IMG]",
    "{H}[IMG, 2]",
    "{H}[IMG, 0.5]",
    "{H}[IMG, {{1.}}]",
    "{H}[IMG, {32, 32}]",
    "{H}[IMG, 2, 0.05]",
    "{H}[IMG, {{0., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}}]",
    '{H}[IMG, "Grayscale"]',
    "{H}[IMG, {1, 0}]",
    "{H}[RAW]",
    # Rank-3 shapes, so "planar only" below means the head really has no volumetric path rather than
    # that no shape fitted it: a 3-D kernel and a 3-element extent.
    "{H}[IMG, {{{1.}}}]",
    "{H}[IMG, {12, 12, 12}]",
]

# BOTH RANKS. The gate would be half a gate covering only planes: two of the three bugs it exists
# for were a volumetric path failing to gain what its planar twin already had, so every head is asked
# the question at rank 2 and at rank 3. 64x64 and 24x24x24 are comfortably above the packing
# threshold and small enough that a few hundred calls finish quickly.
RAW2 = "Table[N[Mod[i*13 + j*7, 97]]/97, {i, 1, 64}, {j, 1, 64}]"
RAW3 = "Table[N[Mod[z*13 + y*7 + x*3, 97]]/97, {z, 1, 24}, {y, 1, 24}, {x, 1, 24}]"
INPUTS = [
    ("2d", 'Image[%s, "Real"]' % RAW2, RAW2),
    ("3d", 'Image3D[%s, "Real"]' % RAW3, RAW3),
]


def heads():
    found = set()
    for rel in SOURCES:
        text = (ROOT / rel).read_text()
        found |= set(re.findall(r'symtab_add_builtin\("([A-Za-z0-9]+)"', text))
    return sorted(found - NOT_IMAGE_RETURNING)


def run(script):
    """Return (stdout, error) -- the error is a string when the probe did not run properly.

    The exit status is CHECKED, and that is not incidental. The first version of this tool returned
    only stdout, and a missing bracket in the generated probe made Mathilda exit 1 with a syntax
    error before printing anything -- whereupon the gate reported "0 image-returning heads" and
    "no newly nested image heads" and exited 0. A gate that passes because its probe never ran is
    worse than no gate: it is a green light with nothing behind it."""
    try:
        out = subprocess.run(
            [str(BIN), "-file", script], capture_output=True, text=True, timeout=120
        )
    except subprocess.TimeoutExpired:
        return "", "the probe timed out"
    err = ""
    if out.returncode != 0:
        err = "the probe exited %d: %s" % (out.returncode,
                                           (out.stderr or out.stdout).strip().splitlines()[-1:])
    return out.stdout, err


def probe(names):
    """One process; one Print per (head, rank, shape): the storage head, or "no"."""
    lines = []
    for h in names:
        for rank, img, raw in INPUTS:
            for i, shape in enumerate(SHAPES):
                call = shape.replace("{H}", h).replace("IMG", img).replace("RAW", raw)
                # No Quiet[]: it is registered but UNIMPLEMENTED, so it stays in the expression and
                # the result arrives as Quiet[NDArray] rather than NDArray. A head that declines
                # simply returns unevaluated, which the "no" branch already covers.
                lines.append(
                    'Print["%s|%s|%d|", Module[{r = %s},'
                    ' If[ImageQ[r] || Image3DQ[r], Head[Part[r, 1]], "no"]]];'
                    % (h, rank, i, call)
                )
    script = ROOT / "tools" / ".image_packing_probe.m"
    script.write_text("\n".join(lines) + "\n")
    out, err = run(str(script))
    if err:
        script.rename(script.with_suffix(".m.failed"))
        print("PROBE FAILED: %s" % err)
        print("the generated probe is kept at %s" % script.with_suffix(".m.failed"))
        return None
    script.unlink(missing_ok=True)

    storage = {}
    for line in out.splitlines():
        m = re.match(r"^([A-Za-z0-9]+)\|(2d|3d)\|(\d+)\|\s*(\S+)\s*$", line.strip())
        if not m:
            continue
        head, rank, res = m.group(1), m.group(2), m.group(4)
        if res in ("NDArray", "List") and (head, rank) not in storage:
            storage[(head, rank)] = res
    return storage


def main():
    if not BIN.exists():
        print("Mathilda is not built; run make first")
        return 0

    names = heads()
    storage = probe(names)
    if storage is None:
        return 1
    if not storage:
        # Every head declining every shape means the shapes are wrong or the binary is broken --
        # either way it is a failure to report, not an empty success.
        print("no head returned an image for any shape; the probe is not testing anything")
        return 1

    packed = sorted(k for k, v in storage.items() if v == "NDArray")
    nested = sorted(k for k, v in storage.items() if v == "List")
    covered = {h for (h, _) in storage}
    no_image = [h for h in names if h not in covered]

    print("%d (head, rank) pairs return an image: %d packed, %d nested"
          % (len(storage), len(packed), len(nested)))
    ranks = {}
    for (h, r) in storage:
        ranks.setdefault(h, []).append(r)
    only2d = sorted(h for h, rs in ranks.items() if rs == ["2d"])
    if only2d:
        print("\nplanar only (no 3-D shape accepted) -- not a failure, but this is where a\n"
              "volumetric path would be missing entirely rather than merely unpacked:")
        print("  " + ", ".join(only2d))
    if no_image:
        print("\nno call shape returned an image (add a shape, or list under NOT_IMAGE_RETURNING):")
        print("  " + ", ".join(no_image))

    new = [k for k in nested if k not in KNOWN_UNPACKED]
    fixed = sorted(k for k in KNOWN_UNPACKED if storage.get(k) == "NDArray")

    if fixed:
        print("\nFIXED -- delete these from KNOWN_UNPACKED:")
        for h, r in fixed:
            print("  %s (%s)" % (h, r))
    if new:
        print("\nNEWLY NESTED -- these return Expr trees where a buffer is expected:")
        for h, r in new:
            print("  %s (%s)" % (h, r))
        print("\nAn image head that returns nested Lists pays for every pixel twice: once to build\n"
              "the nodes and once for whatever reads them. Three such cases have each cost between\n"
              "4x and 23x, with correct answers throughout.")
        return 1
    print("\nno newly nested image heads")
    return 0


if __name__ == "__main__":
    sys.exit(main())
