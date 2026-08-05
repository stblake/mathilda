#!/usr/bin/env python3
"""Tutorial example verifier / transcript generator.

Drives the compiled ``./Mathilda`` exactly like ``generate.py`` does, so the
worked transcripts in ``site/docs/tutorials/*.md`` can be checked against the
real binary (or generated from a list of inputs).

Usage:
  # Run a newline-separated list of input expressions and print In/Out pairs.
  python3 site/verify_tutorial.py transcript inputs.txt

  # Check every ```mathematica``` block in a tutorial: parse its In[k]:= lines,
  # re-run them, and report any Out[k]= line that disagrees with the binary.
  python3 site/verify_tutorial.py check site/docs/tutorials/04-arithmetic.md

A tutorial is run as ONE continuous session (state carries across blocks), and
inputs are matched to outputs by position; the displayed In[k]/Out[k] numbers
are cosmetic and ignored.
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MATHILDA = ROOT / "Mathilda"

OUT_RE = re.compile(r"^Out\[\d+\]=\s?(.*)$")
IN_RE = re.compile(r"^In\[\d+\]:=\s?(.*)$")
OUTLINE_RE = re.compile(r"^Out\[(\d+)\]=\s?(.*)$")
FENCE_RE = re.compile(r"^```")


def run_session(lines, timeout=180):
    """Evaluate `lines` in one Mathilda session and return one output string per
    line, in order.

    Mathilda serves a non-tty stdin with a line-based NDJSON protocol (see
    src/repl.c pipe_mode_loop), NOT the interactive In[]/Out[] transcript — so
    each expression is sent as {"id": k, "expr": ...} and its result read back
    from the matching {"id": k, "type": "expr", "payload": ...} line. The
    payload is exactly the OutputForm string the REPL would print as Out[k]=.
    Matching by id keeps a `;`-suppressed line (payload "Null") and any raw
    Print output from confusing the alignment."""
    import json
    proc = subprocess.Popen([str(MATHILDA)], stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                            text=True, bufsize=1)
    payloads = {}
    try:
        for k, e in enumerate(lines, 1):
            proc.stdin.write(json.dumps({"id": k, "expr": e}) + "\n")
        proc.stdin.write('{"type":"quit"}\n')
        proc.stdin.flush()
        proc.stdin.close()
        for line in proc.stdout:
            line = line.strip()
            if not (line.startswith("{") and '"type"' in line):
                continue  # raw Print/CompilePrint output, or a blank line
            try:
                msg = json.loads(line)
            except ValueError:
                continue
            if msg.get("type") == "expr" and "id" in msg:
                payloads.setdefault(msg["id"], msg.get("payload", ""))
            elif msg.get("type") == "error" and "id" in msg:
                payloads.setdefault(msg["id"], msg.get("message", ""))
        proc.wait(timeout=timeout)
    finally:
        if proc.poll() is None:
            proc.kill()
    # One entry per input line, in order; missing (e.g. a bare-Null with no
    # emit) becomes "" so the position still lines up with `lines`.
    return [payloads.get(k, "") for k in range(1, len(lines) + 1)]


def extract_inputs(md_text):
    """Return list of (input_expr, expected_out) from fenced mathematica blocks.
    expected_out is the single-line Out text following each In (or None)."""
    pairs = []
    in_fence = False
    pending = None
    for line in md_text.splitlines():
        if FENCE_RE.match(line):
            in_fence = not in_fence
            continue
        if not in_fence:
            continue
        mi = IN_RE.match(line)
        mo = OUTLINE_RE.match(line)
        if mi:
            if pending is not None:
                pairs.append((pending, None))
            pending = mi.group(1).strip()
        elif mo and pending is not None:
            pairs.append((pending, mo.group(2).strip()))
            pending = None
    if pending is not None:
        pairs.append((pending, None))
    return pairs


def cmd_transcript(path):
    lines = [l.rstrip("\n") for l in Path(path).read_text().splitlines()
             if l.strip() and not l.lstrip().startswith("#")]
    outs = run_session(lines)
    if len(outs) != len(lines):
        print(f"!! input/output count mismatch: {len(lines)} in, {len(outs)} out",
              file=sys.stderr)
    for i, (expr, out) in enumerate(zip(lines, outs), 1):
        print(f"In[{i}]:= {expr}")
        print(f"Out[{i}]= {out}\n")


def cmd_check(path):
    text = Path(path).read_text()
    pairs = extract_inputs(text)
    inputs = [p[0] for p in pairs]
    outs = run_session(inputs)
    if len(outs) != len(inputs):
        print(f"!! count mismatch: {len(inputs)} inputs, {len(outs)} outputs")
    nbad = 0
    for (expr, expected), actual in zip(pairs, outs):
        actual = actual.strip()
        if expected is None:
            continue  # input with no shown Out[] line — nothing to check
        if expected != actual:
            nbad += 1
            print(f"MISMATCH for: {expr}")
            print(f"   doc:    {expected!r}")
            print(f"   binary: {actual!r}\n")
    print(f"{'OK' if nbad == 0 else 'FAIL'}: {len(inputs)} inputs, {nbad} mismatch(es) in {Path(path).name}")
    return nbad


if __name__ == "__main__":
    if len(sys.argv) < 3 or sys.argv[1] not in ("transcript", "check"):
        print(__doc__)
        sys.exit(2)
    if not MATHILDA.exists():
        sys.exit(f"error: {MATHILDA} not found — run `make` first.")
    if sys.argv[1] == "transcript":
        cmd_transcript(sys.argv[2])
    else:
        sys.exit(1 if cmd_check(sys.argv[2]) else 0)
