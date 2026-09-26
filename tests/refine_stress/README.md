# Refine stress suite

A first-principles regression corpus for `Refine[expr, assum]` and the assumption
engine it shares with `Simplify`, `PossibleZeroQ`, `Element`, and `Assuming`.

- **`corpus.py`** — 159 cases across 12 categories (radicals/powers, Abs/Sign/
  Conj/Arg, Log/Exp, integer trig, Floor/Mod, Re/Im, Element/domains,
  Equal/Unequal, inequalities/CAD, deep positivity, assumption plumbing, and
  adversarial soundness/robustness). Each carries the mathematically-correct
  expected form and a one-line first-principles justification.
- **`run_refine_stress.py`** — the runner.

## Running

```bash
make                                              # build ./Mathilda first
make check-refine-stress                          # or:
python3 tests/refine_stress/run_refine_stress.py            # summary + ratchet
python3 tests/refine_stress/run_refine_stress.py --verbose  # full matrix + every non-PASS case
python3 tests/refine_stress/run_refine_stress.py --timeout 60
```

The runner locates `./Mathilda` relative to itself, so it works from any cwd. It
needs no external `timeout` binary (it uses Python's own per-process timeout) and
no third-party packages.

## How it decides

Each case runs in its **own** `Mathilda -file` process (a hang or crash is
contained, not fatal) and is judged by structural equality **inside** Mathilda:

```
RefineStressAcc = <expr>;
Print["@@PASS@@", TrueQ[RefineStressAcc === (<expected>)]];
Print["@@UNCH@@", TrueQ[RefineStressAcc === (<input>)]];
```

so printer-form differences never cause a false failure. Verdicts: `PASS`,
`UNCHANGED` (result === input — a missing reduction, never a wrong answer),
`DIVERGENT` (neither — a potential wrong answer), `ERROR`, `HANG`, `CRASH`.

## Ratchet

Unlike `tests/test_refine.c` (which asserts the PASSing subset in C), this suite
also holds the cases Refine cannot yet reduce, listed in `BASELINE_UNCHANGED`.
It **fails** when:

- any case is `DIVERGENT` / `ERROR` / `HANG` / `CRASH` (Refine must never return a
  wrong answer, hang, or crash);
- a case **not** on the baseline goes `UNCHANGED` (a reduction regressed);
- a case **on** the baseline now `PASS`es (delete it from the baseline — the
  runner names which);
- the baseline names an id no longer in the corpus.

So improving `Refine` (fixing a baseline gap) turns the suite red until the
baseline line is removed — the same forward-ratchet the `make check-*` audits
use. The standing gaps are documented in `tasks/refine_stress_report.md`.
