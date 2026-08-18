# 87 — Diophantine solving: Mathilda vs sympy

A head-to-head of **integer equation solving** — Mathilda's
`Solve[eqns && constraints, vars, Integers]` (the `src/solveint.c` pre-pass)
against **sympy** (`sympy.solvers.diophantine.diophantine`) — across the
well-known Diophantine families: sums of squares, linear systems, Pell and
negative Pell, Pythagorean triples with a perimeter, Egyptian fractions, Markov
triples, sums of three cubes, Mordell curves, taxicab numbers, Euler's
sum-of-like-powers, Catalan/exponential, and mixed bounded boxes.

## What is measured

For each family the same problem — the *same box* — is put to three solvers:

| Column | What it is |
|--------|------------|
| **Mathilda** | `Solve[..., Integers]`, timed in-process with `AbsoluteTiming` (min of the case's reps). Returns the finite, constraint-satisfying set. |
| **sympy `diophantine`** | The library's actual capability. Classified per run as `solved` (finite comparable set), `parametric` (an unbounded family — no bounded answer), `empty` (ruled out, e.g. by gcd), or `NotImplementedError`. |
| **Python search** | The fallback a sympy user writes when `diophantine` cannot answer: a reasonable search over the *same* box (innermost variable solved in closed form where possible), under a wall-clock budget. Quantifies the efficiency gap. |

Solution **counts** are the invariant compared (ordering and radical spelling
are not). Sum-of-squares cases pass `permute=True` so sympy expands its base
representatives to the full signed/ordered set Mathilda enumerates.

## The three categories

- **A — sympy returns a finite set** (sum of two/three squares). Directly
  comparable; the counts agree exactly (12, 28, 72).
- **B — sympy solves the family but only PARAMETRICALLY** (linear, Pell,
  negative Pell, and — after a hand elimination — the Pythagorean system).
  `diophantine` hands back an infinite parametric family or a raw conic; the
  bounded answer the user asked for still needs a search or filter on top.
- **C — sympy has no solver** (`NotImplementedError`): every mixed/cubic/
  exponential form (motivating `x²+2y³`, Mordell, taxicab, Euler fifths,
  Catalan, boxed powers, Egyptian, Markov, sums of three cubes) **and every
  system** — `diophantine` takes a single equation, so the two-quadratic and
  Pythagorean systems cannot even be expressed.

## Headline

Mathilda answers **all 19** cases with the exact finite solution set (all 19
counts verified). sympy's `diophantine` answers **3** directly and comparably,
handles a handful more only as an unbounded parametric family or after manual
elimination, and raises `NotImplementedError` on **11** — the entire cubic /
exponential / system half of the table. Where a Python search is the only
sympy-ecosystem option, Mathilda is typically orders of magnitude faster and,
for the `10⁵`-box sum of three cubes and the taxicab search, is the only method
that finishes at all.

The current numbers are in [`REPORT.md`](REPORT.md); the machine-readable form
is [`results.json`](results.json).

## Running

```bash
# one-time: a venv with sympy
python3 -m venv ../.venv
../.venv/bin/pip install sympy

# build the Mathilda binary at the repo root first (make -j)
cd benchmarks/87-diophantine-integers
../.venv/bin/python run.py
```

`cases.py` is the single source of truth (one row per family). `run.py`
generates [`diophantine.m`](diophantine.m) from it, runs it once, runs the
sympy and Python columns in-process, joins by label, and writes `REPORT.md` +
`results.json`. Override the binary with `MATHILDA_BIN=...`; the Python-search
budget with `BRUTE_BUDGET_S=...`.

## Held-out validation (the silent-wrong-answer gate)

`cases.py` is *developed-against*: each row was added in the same commit as the
method that solves it, so a green result there measures "the method works on the
example it was written for", not coverage — and cannot see the one failure that
matters, a **silent wrong answer** (a `{}` or finite set that is actually
wrong). [`heldout.py`](heldout.py) is the opposite discipline: ~20 equations
drawn from standard references (Mordell, Nagell, classic lists), **none of them
in `cases.py`**, run cold and cross-checked against an independent Python
brute-force oracle over the same box.

```bash
make check-diophantine-heldout          # or: python3 validate.py
```

[`validate.py`](validate.py) needs only the Mathilda binary (no sympy). For each
case Mathilda emits its status (`uneval` / `empty` / `finite` / `param`) and, for
a concrete answer, the integer solutions inside the box (a parametric family is
expanded over `C[1]`). The verdict is **OK** (matches the oracle), **DECLINE**
(unevaluated where that is acceptable — an honest gap), or **WRONG** (a
`{}`/finite/parametric answer the oracle contradicts). **Any WRONG fails the gate
(nonzero exit).** Results land in [`HELDOUT_REPORT.md`](HELDOUT_REPORT.md). This
gate is what first surfaced the underdetermined-linear-system silent `{}`, and it
now guards every Solve/Integers method against the next such regression.
