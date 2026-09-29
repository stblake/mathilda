#!/usr/bin/env python3
"""check_pipe_protocol.py -- behavioural checks for the NDJSON pipe protocol.

Drives the built ./Mathilda over its non-tty stdin protocol (src/repl.c,
pipe_mode_loop) and asserts on what comes back, in BOTH request modes:

  * plain requests ({"id", "expr"}) -- what the site generator, the audit tools
    and the book tools send. Their behaviour must not move: one expression per
    request, Print text raw on stdout, a Null result sent as "Null".
  * notebook cells ({"id", "expr", "cell": true}) -- what the notebook front
    end sends: several statements per cell, Print output and messages as
    "stream" / "message" lines before the result, no result for `;` or Null.

Plus the request-reader fixes shared by both: a request longer than the old
10 KB line buffer, and JSON \\uXXXX escapes.

    make check-pipe-protocol          (or: python3 tools/check_pipe_protocol.py)

Exit status 1 on any failure; SKIP (exit 0) when ./Mathilda is not built.
"""

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MATHILDA = ROOT / "Mathilda"

failures = 0


def check(name, cond, got=None):
    global failures
    print(f"{'PASS' if cond else 'FAIL'}  {name}" + ("" if cond else f"\n      got: {got!r}"))
    if not cond:
        failures += 1


def session(requests):
    """Send `requests` (dicts) in one session; return (json_by_id, raw_lines, stderr)."""
    lines = [json.dumps(r) for r in requests] + ['{"type":"quit"}']
    proc = subprocess.run([str(MATHILDA)], input="\n".join(lines) + "\n",
                          capture_output=True, text=True, timeout=120,
                          env={"MATHILDA_NO_WINDOW": "1", "PATH": "/usr/bin:/bin"})
    by_id, raw = {}, []
    for line in proc.stdout.splitlines():
        try:
            msg = json.loads(line)
        except ValueError:
            raw.append(line)
            continue
        if isinstance(msg, dict) and "id" in msg:
            by_id.setdefault(msg["id"], []).append(msg)
        elif not (isinstance(msg, dict) and msg.get("type") == "pong"):
            raw.append(line)
    return by_id, raw, proc.stderr


def kinds(msgs):
    return [m["type"] for m in msgs]


def payloads(msgs):
    return [m["payload"] for m in msgs if m["type"] == "expr"]


def main():
    if not MATHILDA.exists():
        print("SKIP  ./Mathilda is not built")
        return 0

    # ---- plain requests: unchanged behaviour -------------------------------
    by_id, raw, err = session([
        {"id": 1, "expr": "1+1"},
        {"id": 2, "expr": "x = 5;"},
        {"id": 3, "expr": "Print[\"@marker\"]; 7"},
        {"id": 4, "expr": "a = 1\nb = 2"},
        {"id": 5, "expr": "1/0"},
    ])
    check("plain: 1+1 -> expr 2 then done", kinds(by_id.get(1, [])) == ["expr", "done"]
          and payloads(by_id[1]) == ["2"], by_id.get(1))
    check("plain: `x = 5;` still sends payload \"Null\" (site/generate.py relies on it)",
          payloads(by_id.get(2, [])) == ["Null"], by_id.get(2))
    check("plain: Print text stays raw on stdout (audit tools read it)",
          "@marker" in raw, raw)
    check("plain: Print sends no stream message", "stream" not in kinds(by_id.get(3, [])),
          by_id.get(3))
    check("plain: one expression per request (multi-line is a parse error)",
          kinds(by_id.get(4, [])) == ["error", "done"], by_id.get(4))
    check("plain: messages stay on stderr", "Power::infy" in err
          and "message" not in kinds(by_id.get(5, [])), (err, by_id.get(5)))

    # ---- notebook cells -----------------------------------------------------
    by_id, raw, err = session([
        {"id": 1, "expr": "a = 1\nb = 2\na + b", "cell": True},
        {"id": 2, "expr": "Print[\"hello\"]; Print[\"world\"]; 1/0", "cell": True},
        {"id": 3, "expr": "x = 5;", "cell": True},
        {"id": 4, "expr": "x = 5; y = 6", "cell": True},
        {"id": 5, "expr": "f[x_] :=\n  x^2\nf[3]", "cell": True},
        {"id": 6, "expr": "1 + 1\nf[1,", "cell": True},
        {"id": 7, "expr": "(* just a comment *)", "cell": True},
        {"id": 8, "expr": "Print[1]\nSin[1, 2]\n3", "cell": True},
        {"id": 9, "expr": "?Sin", "cell": True},
    ])
    check("cell: one statement per line, every result shown",
          payloads(by_id.get(1, [])) == ["1", "2", "3"], by_id.get(1))
    m2 = by_id.get(2, [])
    check("cell: Print output arrives as stream lines, before the result",
          kinds(m2) == ["stream", "stream", "message", "expr", "done"]
          and "".join(m["text"] for m in m2 if m["type"] == "stream") == "hello\nworld\n", m2)
    check("cell: a message arrives as a message line",
          any(m["type"] == "message" and m["text"].startswith("Power::infy") for m in m2), m2)
    check("cell: nothing raw on stdout", raw == [], raw)
    check("cell: nothing on stderr during evaluation", "Power::infy" not in err, err)
    check("cell: `x = 5;` sends no result", kinds(by_id.get(3, [])) == ["done"], by_id.get(3))
    check("cell: `x = 5; y = 6` shows only y's value",
          payloads(by_id.get(4, [])) == ["6"], by_id.get(4))
    check("cell: a statement continues over a line break inside a definition",
          payloads(by_id.get(5, [])) == ["9"], by_id.get(5))
    check("cell: a syntax error anywhere evaluates nothing",
          kinds(by_id.get(6, [])) == ["error", "done"], by_id.get(6))
    check("cell: a comment-only cell is just done", kinds(by_id.get(7, [])) == ["done"],
          by_id.get(7))
    m8 = by_id.get(8, [])
    check("cell: per-statement order (Print, then the next statement's message and result)",
          kinds(m8) == ["stream", "message", "expr", "expr", "done"], m8)
    check("cell: ?Sin still answers as a usage message",
          kinds(by_id.get(9, [])) == ["usage", "done"], by_id.get(9))

    # ---- request reader: length and escapes ---------------------------------
    big = "Length[{" + ",".join(["1"] * 20000) + "}]"          # ~40 KB, past the old 10 KB buffer
    by_id, raw, err = session([
        {"id": 1, "expr": big},
        {"id": 2, "expr": big, "cell": True},
        {"id": 3, "expr": "StringLength[\"a\tb\u0001c\"]"},   # serde writes \u0001
        {"id": 4, "expr": "\"é\U0001F600\""},                # json.dumps escapes both
    ])
    check("a 40 KB plain request is read whole", payloads(by_id.get(1, [])) == ["20000"],
          by_id.get(1))
    check("a 40 KB cell is read whole", payloads(by_id.get(2, [])) == ["20000"], by_id.get(2))
    check("\\u00XX control escapes decode", payloads(by_id.get(3, [])) == ["5"], by_id.get(3))
    check("\\uXXXX and surrogate pairs decode to UTF-8",
          payloads(by_id.get(4, [])) == ["\"é\U0001F600\""], by_id.get(4))

    print(f"\n{'OK' if not failures else f'{failures} FAILED'}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
