# Refine stress-test report (v0.202 baseline → v0.203)

A first-principles corpus of **159 `Refine` cases** across 12 categories was built to test every
documented surface of `Refine[expr, assum]` plus adversarial soundness/robustness probes. Each
expected answer is derived from mathematical first principles (no Mathematica oracle); the harness
categorizes by structural equality (`SameQ`) so printer-form differences never cause a false
failure. Corpus + driver live in the session scratchpad
(`refine_stress/{corpus.py,run_stress.py}`); each case runs in its own process under a hard timeout
so hangs/crashes are first-class outcomes.

## Verdicts

| Verdict | Meaning |
|---|---|
| PASS | result `===` the first-principles expected value |
| UNCHANGED | result `===` the input — a missing reduction (feature gap), never a wrong answer |
| DIVERGENT | neither — hand-reviewed for soundness vs sound-but-incomplete |
| HANG / CRASH | process timeout / died by signal |

## Result: 137/159 → 151/159, zero wrong answers

| category | baseline PASS | after v0.203 |
|---|---|---|
| 1 Radicals/Powers | 15 | **18** |
| 2 Abs/Sign/Conj/Arg | 14 | **17** |
| 3 Log/Exp | 6 | **9** |
| 4 Integer trig/parity | 9 | 9 |
| 5 Floor/Ceil/Mod | 13 | 13 |
| 6 Re/Im/ComplexExpand | 9 | **10** |
| 7 Element/domain | 11 | **15** |
| 8 Equal/Unequal | 9 | **10** |
| 9 Inequalities/logic | 15 | 15 |
| 10 Deep positivity | 6 | 6 |
| 11 Plumbing | 14 | 14 |
| 12 Adversarial/soundness | 14 | **15** |
| **TOTAL** | **137** | **151** |

Final: **151 PASS / 8 UNCHANGED / 0 DIVERGENT / 0 ERROR / 0 HANG / 0 CRASH.** Every remaining
non-pass returns the input unchanged (a missing reduction) — there are no incorrect answers.

## What the corpus caught

### Soundness bug (fixed) — highest priority
- **`Refine[a == b, a - b == 0]` → `False`** (should be `True`), and even
  **`Refine[a == b, a == b]` → `False`**; **`PossibleZeroQ[a - b, Assumptions -> a - b == 0]` →
  `False`**. Root cause: the assumption-aware Schwartz–Zippel sampler (`src/zero_test.c`,
  `extract_spec`) folds only single-variable constraints and ignores a *coupling* equality, then
  samples points the assumptions exclude and reports a genuine identity as non-zero. **Fix:**
  downgrade a sampler `False` to `Unknown` when a ctx fact touching the expression's symbols is an
  equality/compound membership it cannot honour (inequalities keep `False` sound — full-measure
  region); and in `refine.c`'s `Equal` branch, prove `True` by equality-substitution. Repairs both
  `Refine` and `PossibleZeroQ`.

### Missing-feature gaps (fixed)
| ids | gap | fix |
|---|---|---|
| rad05, abs07 | `x <= 0` (NonPositive) not reduced (`Sqrt[x^2]`/`Abs` gave `Abs[x]`) | NonPositive rule bucket in `simp_assume_rewrite.c` |
| abs10, abs11 | `Arg[x]` under sign facts | `Arg[x] -> 0 / Pi` string rules |
| cpx07 | `Abs[a+bI]` real parts → magnitude | add `Abs` to the ComplexExpand-under-reality pass |
| elt10–12, elt16 | `Element[x, Positive/Negative/NonNegative]` as queried domain | delegate to sign provers in `element_decide` |
| log03, log05, log06 | `Log[x^2]→2Log|x|`, `Log[E^x]→x` (real), `Log[x·rest]→Log[x]+Log[rest]` (pos) | string rules |
| adv06 | 16+ assumed symbols overflow the rule buffer → **all** rules dropped | prune rule synthesis to symbols in the target |

### Methodology check — the corpus caught *my* mistakes
Two initial expected values were wrong and were corrected (Mathilda was right): `plm03`
(`-Sqrt[x^2] y`: with only `y<0`, `x` is not even real, so `√(x²)` must stay) and `dpp02`
(`Sign[x^2+1]` with no assumption stays — `x` could be complex). One harness bug was also caught:
the accumulator variable `r` collided with the data symbol `r` in `(x^2)^r`, self-referencing into
the recursion limit (a false HANG) — a known trap; the accumulator was renamed.

## Open follow-ups (return input unchanged, never wrong — not attempted here)

| ids | gap | why deferred |
|---|---|---|
| trg08, trg09, trg10 | integer-*linear* trig args: `Exp[2 Pi I k]→1`, `Cos[2 k Pi]→1`, `Sin[(2k+1)Pi/2]→(-1)^k` | rules key on the bare integer symbol `k`; needs coefficient-aware matching |
| elt09 | `Element[Sqrt[2], Algebraics] → True` | needs recognition of algebraic literals (radicals/`Root`/`AlgebraicNumber`) |
| dpp03, dpp04 | `Abs[(x-1)^2]→(x-1)^2`, `Sqrt[(x-1)^2]→Abs[x-1]` | deep-positivity pass proves only *strict* signs and uses domain Reals; extending to non-strict / compound-real bases needs a reality check it doesn't currently carry |
| inq16 | 7-variable entailment | CAD bails at `> 6` bare variables (`reduce_is_unsat`) |
| adv03 | `Sign[(√2+√3)^2 - 5 - 2√6] → 0` | argument is an exact algebraic zero; needs zero-normalization of the argument before `Sign` (currently sound: leaves it, never ±1) |

## How to re-run

```bash
make -j$(nproc)                       # build ./Mathilda
python3 <scratch>/refine_stress/run_stress.py 40   # full corpus, per-case timeout 40s
cd tests/build && cmake .. && make refine_tests && ./refine_tests
```
