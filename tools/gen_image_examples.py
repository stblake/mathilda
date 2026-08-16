#!/usr/bin/env python3
"""gen_image_examples.py -- write topic-grouped examples for every image head.

WHY A GENERATOR. Every image page needs the same subsections (Basic Examples, Scope, Options,
Applications, Properties & Relations, Neat Examples) and enough examples in each to be worth
reading. Hand-writing that for thirty heads is thirty chances to write an example that does not
run, and the failure is silent: `site/generate.py` re-verifies every example against the binary,
so a wrong one simply loses its output and the page shows a call with no answer.

So candidates are PROPOSED here and only the ones that actually evaluate are written out. The
proposal is per-family rather than per-head -- a filter taking a radius accepts the same shapes
whichever filter it is -- and the verification is per-candidate, which is what keeps a family
template from putting a call on a page that does not support it.

WHAT MAKES AN EXAMPLE GOOD. Variety, not volume: the same filter at six radii on one image is
six rows that say one thing. Each group draws on a bank of distinctly different images (a
checkerboard, a disk, a lighting ramp, a zone plate, a colour gradient, a bit image, a byte
image, a volume) so that the reader sees what a head does to structure, to smooth gradients and
to noise -- and the identities in Properties & Relations state the algebra that a picture cannot
show.

EACH FENCE IS ITS OWN SESSION. `site/generate.py` verifies one fenced block per session, so a
block that uses a name must define it. Every emitted block therefore opens with its own setup
line. Getting this wrong is how six examples on the Import page came back unevaluated.
"""

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "Mathilda"
SPEC = ROOT / "docs" / "spec" / "builtins" / "image-processing.md"

# ---------------------------------------------------------------- the image bank
#
# Each entry is (name, definition, one-line description). The description becomes the example's
# prose, so a reader knows WHICH image is being shown rather than seeing eight anonymous calls.
IMAGES = {
    "chk": (
        'chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];',
        "a checkerboard: hard edges in both directions",
    ),
    "disk": (
        'disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];',
        "a disk: one closed curved boundary",
    ),
    "ramp": (
        'ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];',
        "a linear ramp: no structure, only a gradient",
    ),
    "zone": (
        'zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];',
        "a zone plate: every spatial frequency in one picture",
    ),
    "noise": (
        'noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];',
        "a deterministic pseudo-noise field",
    ),
    "rgb": (
        'rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];',
        "a three-channel colour gradient",
    ),
    "sky": (
        'sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];',
        "a sky-like vertical colour gradient",
    ),
    "bit": (
        "bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];",
        'a "Bit" image, whose stored values are 0 and 1',
    ),
    "byte": (
        "byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];",
        'a "Byte" image, whose stored values run 0..255',
    ),
    "vol": (
        'vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];',
        "a volume, 12 x 10 x 8",
    ),
    "volb": (
        'volb = Image3D[Table[N[Boole[x <= 6 && y <= 5]], {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];',
        "a volume with one straight interior edge",
    ),
}

# ------------------------------------------------------------------- families
#
# A family is a set of call shapes plus the identities that hold for it. `{H}` is the head and
# `{v}` an image variable.
FILTER_LIKE = {
    "basic": [
        "{H}[chk]",
        "{H}[chk, 2]",
        "ImageDimensions[{H}[chk, 2]]",
        "ImageType[{H}[chk, 1]]",
    ],
    "scope": [
        "{H}[disk, 1]",
        "{H}[ramp, 2]",
        "{H}[zone, 2]",
        "{H}[noise, 3]",
        "{H}[rgb, 1]",
        "{H}[sky, 2]",
        "{H}[bit, 1]",
        "{H}[byte, 2]",
        "{H}[vol, 1]",
        "ImageChannels[{H}[rgb, 2]]",
        "ImageDimensions[{H}[vol, 1]]",
        "{H}[chk, {1, 3}]",
        "{H}[chk, 4]",
    ],
    "rel": [
        "ImageDimensions[{H}[chk, 3]] === ImageDimensions[chk]",
        "ImageChannels[{H}[rgb, 2]] === ImageChannels[rgb]",
        "ImageData[{H}[ramp, 0]] === ImageData[ramp]",
        "Max[Flatten[ImageData[{H}[chk, 2]]]] <= 1.0",
        "Min[Flatten[ImageData[{H}[chk, 2]]]] >= 0.0",
        "ImageDimensions[{H}[vol, 2]] === ImageDimensions[vol]",
    ],
    "app": [
        "Binarize[{H}[noise, 2]]",
        "EdgeDetect[{H}[zone, 2]]",
        'ImageDimensions[{H}[Import[Export["/tmp/mathilda_ex.png", rgb]], 2]]',
    ],
    "neat": ["{H}[zone, 4]", "{H}[zone, 1]"],
    "options": [
        "{H}[chk, 2, Padding -> 0]",
        '{H}[chk, 2, Padding -> "Fixed"]',
        '{H}[chk, 2, Padding -> "Reversed"]',
        '{H}[chk, 2, Method -> "Box"]',
    ],
}

MORPH = {
    "basic": [
        "{H}[disk, 1]",
        "{H}[disk, 2]",
        "ImageDimensions[{H}[disk, 2]]",
        "{H}[bit, 1]",
    ],
    "scope": [
        "{H}[chk, 1]",
        "{H}[ramp, 1]",
        "{H}[noise, 2]",
        "{H}[rgb, 1]",
        "{H}[byte, 1]",
        "{H}[vol, 1]",
        "{H}[volb, 1]",
        "{H}[disk, 3]",
        "{H}[disk, 4]",
        "ImageChannels[{H}[rgb, 1]]",
        "ImageDimensions[{H}[vol, 2]]",
    ],
    "rel": [
        "ImageData[{H}[disk, 0]] === ImageData[disk]",
        "ImageDimensions[{H}[disk, 3]] === ImageDimensions[disk]",
        "ImageData[{H}[{H}[disk, 1], 1]] === ImageData[{H}[disk, 2]]",
        "Max[Flatten[ImageData[{H}[bit, 1]]]] <= 1.0",
    ],
    "app": ["Binarize[{H}[noise, 1]]", "EdgeDetect[{H}[disk, 1]]"],
    "neat": ["{H}[zone, 2]"],
    "options": [],
}

GEOM = {
    "basic": ["{H}[chk, 2]", "ImageDimensions[{H}[chk, 2]]", "{H}[disk, 1]"],
    "scope": [
        "{H}[rgb, 2]",
        "{H}[bit, 1]",
        "{H}[byte, 2]",
        "{H}[vol, 1]",
        "{H}[chk, {2, 4}]",
        "{H}[chk, {{1, 2}, {3, 4}}]",
        "ImageDimensions[{H}[vol, 1]]",
        "ImageChannels[{H}[rgb, 2]]",
    ],
    "rel": [
        "ImageChannels[{H}[rgb, 2]] === ImageChannels[rgb]",
        "ImageDimensions[{H}[chk, 0]] === ImageDimensions[chk]",
    ],
    "app": ["EdgeDetect[{H}[disk, 2]]"],
    "neat": ["{H}[zone, 4]"],
    "options": ["{H}[chk, 2, Padding -> 0]", "{H}[chk, 2, Padding -> 1]"],
}

QUERY = {
    "basic": ["{H}[chk]", "{H}[rgb]", "{H}[bit]", "{H}[byte]"],
    "scope": [
        "{H}[disk]",
        "{H}[ramp]",
        "{H}[zone]",
        "{H}[noise]",
        "{H}[sky]",
        "{H}[vol]",
        "{H}[volb]",
        '{H}[Import[Export["/tmp/mathilda_ex.png", rgb]]]',
    ],
    "rel": [
        "{H}[chk] === {H}[GaussianFilter[chk, 1]]",
        "{H}[rgb] === {H}[ImagePad[rgb, 2]]",
    ],
    "app": ["Table[{H}[GaussianFilter[chk, r]], {r, 1, 3}]"],
    "neat": ["{H}[zone]"],
    "options": [],
}

PLAIN = {  # a head taking just an image, with no radius
    "basic": ["{H}[chk]", "{H}[disk]", "ImageDimensions[{H}[chk]]"],
    "scope": [
        "{H}[ramp]",
        "{H}[zone]",
        "{H}[noise]",
        "{H}[rgb]",
        "{H}[sky]",
        "{H}[bit]",
        "{H}[byte]",
        "{H}[vol]",
        "{H}[volb]",
        "ImageChannels[{H}[rgb]]",
        "ImageDimensions[{H}[vol]]",
    ],
    "rel": [
        "ImageDimensions[{H}[chk]] === ImageDimensions[chk]",
        "Max[Flatten[ImageData[{H}[chk]]]] <= 1.0",
        "Min[Flatten[ImageData[{H}[chk]]]] >= 0.0",
    ],
    "app": ["Binarize[{H}[zone]]", "Dilation[{H}[disk], 1]"],
    "neat": ["{H}[zone]"],
    "options": ["{H}[chk, Padding -> 0]", '{H}[chk, Method -> "Box"]'],
}

FAMILY = {
    "GaussianFilter": FILTER_LIKE,
    "MeanFilter": FILTER_LIKE,
    "MedianFilter": FILTER_LIKE,
    "Sharpen": FILTER_LIKE,
    "GradientFilter": FILTER_LIKE,
    "CornerFilter": FILTER_LIKE,
    "LocalAdaptiveBinarize": FILTER_LIKE,
    "Dilation": MORPH,
    "Erosion": MORPH,
    "Opening": MORPH,
    "Closing": MORPH,
    "ImagePad": GEOM,
    "ImageCrop": GEOM,
    "ImageRotate": GEOM,
    "ImageResize": GEOM,
    "ImageQ": QUERY,
    "Image3DQ": QUERY,
    "ImageDimensions": QUERY,
    "ImageChannels": QUERY,
    "ImageType": QUERY,
    "ImageData": QUERY,
    "ImageLevels": QUERY,
    "FindThreshold": QUERY,
    "Binarize": PLAIN,
    "EdgeDetect": PLAIN,
    "ColorNegate": PLAIN,
    "ImageReflect": PLAIN,
    "DistanceTransform": PLAIN,
    "ImageAdjust": PLAIN,
    "ImageCorners": PLAIN,
    "ColorConvert": PLAIN,
    "DerivativeFilter": PLAIN,
}

GROUPS = [
    ("Basic Examples", "basic"),
    ("Scope", "scope"),
    ("Options", "options"),
    ("Applications", "app"),
    ("Properties & Relations", "rel"),
    ("Neat Examples", "neat"),
]

VAR_RE = re.compile(r"\b(chk|disk|ramp|zone|noise|rgb|sky|bit|byte|vol|volb)\b")


def setup_for(exprs):
    """The definitions an expression list needs, in bank order, one per line."""
    used = set()
    for e in exprs:
        used |= set(VAR_RE.findall(e))
    return [IMAGES[k][0] for k in IMAGES if k in used]


def probe(cands):
    """Run every candidate and return {id: printed result}. One process, one Print per candidate.

    Each candidate gets its own Module so a failure cannot leak a partial definition into the
    next, and every image variable is defined once at the top -- this is a probe, not the page,
    so it does not have to obey the one-session-per-fence rule the emitted blocks do."""
    lines = [d for d, _ in IMAGES.values()]
    for cid, expr in cands:
        lines.append('Print["%s|", %s];' % (cid, expr))
    script = ROOT / "tools" / ".image_examples_probe.m"
    script.write_text("\n".join(lines) + "\n")
    r = subprocess.run(
        [str(BIN), "-file", str(script)], capture_output=True, text=True, timeout=900
    )
    script.unlink(missing_ok=True)
    out = {}
    for line in r.stdout.splitlines():
        m = re.match(r"^([A-Za-z0-9]+_[a-z]+_\d+)\|\s*(.*)$", line)
        if m:
            out[m.group(1)] = m.group(2).strip()
    return out


def acceptable(head, expr, res):
    """Is this a result worth putting on a page?"""
    if not res:
        return False
    if res.startswith(head + "["):  # declined: the call came back unevaluated
        return False
    for bad in ("$Failed", "Indeterminate", "ComplexInfinity", "::", "Null"):
        if bad in res:
            return False
    if len(res) > 300:  # a wall of numbers teaches nothing
        return False
    # An identity example is only worth showing when it is TRUE: a False would document a bug
    # rather than a property, and if one appears the property is wrong and should be fixed here.
    if (
        "===" in expr
        or expr.startswith(("Max[", "Min["))
        or " <= " in expr
        or " >= " in expr
    ):
        return res == "True"
    return True


def main():
    heads = sorted(FAMILY)
    cands = []
    for h in heads:
        fam = FAMILY[h]
        for _, key in GROUPS:
            for n, tpl in enumerate(fam[key]):
                cands.append(("%s_%s_%d" % (h, key, n), tpl.replace("{H}", h)))
    print("probing %d candidates over %d heads..." % (len(cands), len(heads)))
    res = probe(cands)
    if not res:
        print("the probe returned nothing; refusing to rewrite the spec")
        return 1

    kept = {}
    for cid, expr in cands:
        h, key, _ = cid.split("_")
        r = res.get(cid, "")
        if acceptable(h, expr, r):
            kept.setdefault((h, key), []).append((expr, r))

    # ---- emit ----------------------------------------------------------------
    text = SPEC.read_text()
    total_new = 0
    for h in heads:
        m = re.search(r"^## %s\s*$" % re.escape(h), text, re.M)
        if not m:
            continue
        # The end of this section: the next H2, or end of file.
        nxt = re.search(r"^## \w+\s*$", text[m.end() :], re.M)
        end = m.end() + (nxt.start() if nxt else len(text) - m.end())

        blocks = []
        for title, key in GROUPS:
            pairs = kept.get((h, key), [])
            if not pairs:
                continue
            setup = setup_for([e for e, _ in pairs])
            body = ["#### %s" % title, "", "```mathematica"]
            n = 1
            for line in setup:
                body.append("In[%d]:= %s" % (n, line))
                body.append("")
                n += 1
            for expr, r in pairs:
                body.append("In[%d]:= %s" % (n, expr))
                body.append("Out[%d]= %s" % (n, r))
                body.append("")
                n += 1
            if body[-1] == "":
                body.pop()
            body.append("```")
            blocks.append("\n".join(body))
            total_new += len(pairs)
        if blocks:
            text = (
                text[:end].rstrip() + "\n\n" + "\n\n".join(blocks) + "\n\n" + text[end:]
            )

    SPEC.write_text(text)
    print("wrote %d verified examples" % total_new)
    per = {}
    for (h, _), v in kept.items():
        per[h] = per.get(h, 0) + len(v)
    for h in sorted(per, key=lambda k: -per[k]):
        print("  %-24s %d" % (h, per[h]))
    dropped = len(cands) - total_new
    print(
        "\n%d candidates dropped (declined, failed, or an identity that came back False)"
        % dropped
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
