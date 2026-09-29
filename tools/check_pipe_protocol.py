#!/usr/bin/env python3
"""check_pipe_protocol.py -- behavioural checks for the NDJSON pipe protocol.

Drives the built ./Mathilda over its non-tty stdin protocol (src/repl.c,
pipe_mode_loop) and asserts on what comes back, in BOTH request modes:

  * plain requests ({"id", "expr"}) -- what the site generator, the audit tools
    and the book tools send. Their behaviour must not move: one expression per
    request, Print text raw on stdout, a Null result sent as "Null".
  * notebook cells ({"id", "expr", "cell": true}) -- what the notebook front
    end sends: several statements per cell, Print output and messages as
    "stream" / "message" lines before the result, no result for `;` or Null,
    and a "line" line per statement carrying the $Line it took, which is what
    makes `%` / `%%` / `%n` / In[n] / Out[n] work and what the front end labels
    the cell In[n] with.

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


def lines(msgs):
    """The $Line numbers a cell reported, in order."""
    return [m["line"] for m in msgs if m["type"] == "line"]


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
    check("plain: no line numbering -- a batch tool has no history, and Out[n] "
          "would hold every result of a sweep alive",
          all("line" not in kinds(m) for m in by_id.values()), by_id)

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
          kinds(m2) == ["line", "stream", "line", "stream", "line", "message", "expr", "done"]
          and "".join(m["text"] for m in m2 if m["type"] == "stream") == "hello\nworld\n", m2)
    check("cell: a message arrives as a message line",
          any(m["type"] == "message" and m["text"].startswith("Power::infy") for m in m2), m2)
    check("cell: nothing raw on stdout", raw == [], raw)
    check("cell: nothing on stderr during evaluation", "Power::infy" not in err, err)
    check("cell: `x = 5;` sends no result, only its line",
          kinds(by_id.get(3, [])) == ["line", "done"], by_id.get(3))
    check("cell: `x = 5; y = 6` shows only y's value",
          payloads(by_id.get(4, [])) == ["6"], by_id.get(4))
    check("cell: a statement continues over a line break inside a definition",
          payloads(by_id.get(5, [])) == ["9"], by_id.get(5))
    check("cell: a syntax error anywhere evaluates nothing",
          kinds(by_id.get(6, [])) == ["error", "done"], by_id.get(6))
    check("cell: a comment-only cell is just done", kinds(by_id.get(7, [])) == ["done"],
          by_id.get(7))
    m8 = by_id.get(8, [])
    check("cell: per-statement order (line, Print, then the next statement's line, message "
          "and result)",
          kinds(m8) == ["line", "stream", "line", "message", "expr", "line", "expr", "done"], m8)
    check("cell: ?Sin still answers as a usage message",
          kinds(by_id.get(9, [])) == ["line", "usage", "done"], by_id.get(9))

    # ---- cell session history: %, %%, %n, In[n], Out[n] ---------------------
    #
    # `%` is not a REPL feature. It parses to Out[-1], which builtin_out resolves
    # against $Line, so a front end that evaluates without recording $Line /
    # In[n] / Out[n] leaves it unevaluated -- which is what the notebook did:
    # `Integrate[x^5 E^x, x]` then `% // Factor` printed a literal `Out[-1]`.
    by_id, raw, err = session([
        {"id": 1, "expr": "11 + 11", "cell": True},
        {"id": 2, "expr": "% + 1", "cell": True},
        {"id": 3, "expr": "%%", "cell": True},
        {"id": 4, "expr": "Out[1]", "cell": True},
        {"id": 5, "expr": "1 + 1\n% * 10", "cell": True},
        {"id": 6, "expr": "z = 99;", "cell": True},
        {"id": 7, "expr": "%", "cell": True},
        {"id": 8, "expr": "$Line", "cell": True},
    ])
    check("cell: `%` is the previous result, not a literal Out[-1]",
          payloads(by_id.get(2, [])) == ["23"], by_id.get(2))
    check("cell: `%%` reaches two results back", payloads(by_id.get(3, [])) == ["22"],
          by_id.get(3))
    check("cell: Out[n] addresses a numbered line", payloads(by_id.get(4, [])) == ["22"],
          by_id.get(4))
    check("cell: `%` inside one cell is the line above, not the cell above",
          payloads(by_id.get(5, [])) == ["2", "20"], by_id.get(5))
    check("cell: a `;`-suppressed statement still counts, and `%` still finds its value",
          payloads(by_id.get(7, [])) == ["99"], by_id.get(7))
    check("cell: one line per STATEMENT, numbered from 1 and never reused",
          [lines(by_id.get(i, [])) for i in range(1, 9)]
          == [[1], [2], [3], [4], [5, 6], [7], [8], [9]],
          [lines(by_id.get(i, [])) for i in range(1, 9)])
    check("cell: $Line is the line being evaluated", payloads(by_id.get(8, [])) == ["9"],
          by_id.get(8))

    # ---- the LaTeX the notebook renders ------------------------------------
    #
    # The front end typesets `latex` when it is present, so a head the LaTeX
    # printer does not know shows as FullForm dressed up as mathematics. Every
    # Solve/DSolve answer is built from Rule, which is how that surfaced.
    def latex(msgs):
        return [m.get("latex") for m in msgs if m["type"] == "expr"]

    by_id, raw, err = session([
        {"id": 1, "expr": "a -> b", "cell": True},
        {"id": 2, "expr": "a :> b", "cell": True},
        {"id": 3, "expr": "Solve[x^2 == 4, x]", "cell": True},
        {"id": 4, "expr": "Reduce[x^2 > 1, x]", "cell": True},
        {"id": 5, "expr": "1 < x < 2", "cell": True},
        {"id": 6, "expr": "!a", "cell": True},
        {"id": 7, "expr": "True", "cell": True},
    ])
    check("latex: Rule prints infix, not Rule[a, b]", latex(by_id.get(1, [])) == ["a\\to b"],
          by_id.get(1))
    check("latex: RuleDelayed is distinguishable from Rule",
          latex(by_id.get(2, [])) == ["a:\\to b"], by_id.get(2))
    check("latex: a Solve answer is readable mathematics",
          latex(by_id.get(3, [])) == ["\\{\\{x\\to -2\\}, \\{x\\to 2\\}\\}"], by_id.get(3))
    check("latex: Or/relations print infix", latex(by_id.get(4, [])) == ["x<-1\\lor x>1"],
          by_id.get(4))
    check("latex: the Inequality chain prints as a chain",
          latex(by_id.get(5, [])) == ["1<x<2"], by_id.get(5))
    check("latex: Not prints as a negation", latex(by_id.get(6, [])) == ["\\neg a"],
          by_id.get(6))
    check("latex: True is a word, not a product of italic letters",
          latex(by_id.get(7, [])) == ["\\text{True}"], by_id.get(7))

    # ---- the printer directives: no latex at all ----------------------------
    #
    # InputForm, FullForm, TeXForm and NumberForm ask for a notation that is not
    # StandardForm, and typesetting has only StandardForm to give. The kernel
    # must therefore send NO `latex` field, which is how the front end knows to
    # show the payload as text. Sending one is not a cosmetic slip: the notebook
    # prefers `latex` whenever it is there, so `expr // InputForm` typeset the
    # StandardForm the reader had just asked not to see, and the directive did
    # nothing outside the terminal REPL. A directive nested inside the result
    # counts too -- the check is over the whole tree, not the outer head.
    by_id, raw, err = session([
        {"id": 1, "expr": "D[Log[1 - Sqrt[x]] Sqrt[x], x] // InputForm", "cell": True},
        {"id": 2, "expr": "FullForm[a + b]", "cell": True},
        {"id": 3, "expr": "TeXForm[a/b]", "cell": True},
        {"id": 4, "expr": "NumberForm[1.23456789, 4]", "cell": True},
        {"id": 5, "expr": "Hold[InputForm[x + y]]", "cell": True},
        {"id": 6, "expr": "D[Log[1 - Sqrt[x]] Sqrt[x], x]", "cell": True},
    ])
    for i, name in ((1, "InputForm"), (2, "FullForm"), (3, "TeXForm"),
                    (4, "NumberForm"), (5, "a nested InputForm")):
        check(f"latex: {name} sends no latex, so the payload is shown as text",
              latex(by_id.get(i, [])) == [None], by_id.get(i))
    check("latex: InputForm's payload is the input form, not StandardForm",
          payloads(by_id.get(1, [])) == ["-1/2/(1 - Sqrt[x]) + 1/2 Log[1 - Sqrt[x]]/Sqrt[x]"],
          by_id.get(1))
    # The same derivative WITHOUT the directive still typesets: the scan for a
    # directive must not cost an ordinary result its LaTeX.
    check("latex: the undirected result still typesets",
          latex(by_id.get(6, [])) not in ([None], []), by_id.get(6))

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
