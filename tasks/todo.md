# Seventh round — closing the coverage sweep's NumPy gaps (2026-08-01)

Working from `plans/HPC_IMPROVEMENT_PLAN.md` §11 (the coverage sweep's register,
experiment 20). Every item below is a **measured** row from
`tools/numeric_sweep.py`, re-confirmed on this machine before starting.

Previous rounds: `docs/design/performance.md` §8–§13, `docs/experiments/`.

Re-measured 2026-08-01, before any change:

| probe | Mathilda | NumPy | ratio |
|---|---:|---:|---:|
| `Extract[v, {{1},{2},{3}}]` | 99.7 ms | 1.0 µs | 97139× |
| `MatrixPower[A300, 4]` | 14.7 s | 865 µs | ~17000× |
| `LeastSquares[A500, b500]` | **does not finish in 180 s** | — | — |
| `PseudoInverse[A300]` | **does not finish** | — | — |
| `Subdivide[0., 1., 999999]` | 1.87 s | 1.08 ms | 1740× |
| `Array[N, 10^6]` | 446 ms | 468 µs | 953× |
| `Table[N[i], {i, 10^6}]` | 428 ms | 474 µs | 902× |
| `UnitVector[10^6, 1]` | 130 ms | 267 µs | 487× |
| `IdentityMatrix[1000]` | 77.3 ms | 231 µs | 334× |
| `DiagonalMatrix[b1000]` | 68.0 ms | 231 µs | 294× |

## The finding that orders the work

`Extract`, `MatrixPower`, `LeastSquares` and `PseudoInverse` are **already on
`pack.c`'s `AWARE` list**, so the transparency gate does *not* materialise for
them. They are slow for a different reason, and it is the same reason in all
four: the builtin itself declines the buffer on its first line —
`linalg_delist_and_reeval(res)`, or (for `Extract`) `expr_part` rejecting an
atomic NDArray and the **post-gate** materialising on the way to rest. The
answer is right, the buffer is thrown away, and a 300×300 machine matrix runs
the exact symbolic path.

`MATHILDA_PACK_DIAG=gate` is what found it — the *post*-gate report names the
head. It is a fifth instance of this document's recurring defect in a new
disguise: not "the head is not aware", but "the head is aware and declines
anyway".

## Items

- [x] **C.11a `LeastSquares` on a machine matrix does not terminate.** Fixed via
      one thin `gesdd` (not `dgels`: the SVD handles rank deficiency, which is
      what `Direct`'s documented "minimum-norm minimiser" requires).
- [x] **C.11b `PseudoInverse` likewise** — same SVD, `A⁺ = V Σ⁺ Uᵀ`.
- [x] **C.2 `MatrixPower`** — binary exponentiation over `dot2`.
- [x] **C.4 `Extract`** — `ndarray_part`, and `INT64_OK`.
- [x] **C.5 the producers** — `Subdivide`, `IdentityMatrix`, `UnitVector`,
      `DiagonalMatrix`, `Rescale` done. **`Array` and `Table` with an exact
      iterator deliberately NOT done** — see Review.
- [x] **Not on the list: `ProductLog`'s machine kernel returned a wrong value**
      for `1 < x < ~1.05`. Found while diagnosing C.3.
- [x] **Not on the list: `ConjugateTranspose` was 195× NumPy.** Turned up while
      re-measuring the linalg group.
- [x] Verify: differential sweep, full test suite, bench gates, valgrind
- [x] Document: `performance.md` §14, weekly changelog, HPC plan register,
      six `docs/spec/builtins/` pages, `tasks/lessons.md`

## Non-negotiables carried from the plan

- **N1** — never trade an exact answer for a fast one. An integer producer
  writes `NDT_INT64`, never float64.
- **N2** — representation is never observable. Every item ships with a
  differential sweep: the same expression packed and unpacked, printed forms
  byte-identical.

## Review

### Results

One uncontended `tools/numeric_sweep.py` run, so both columns are the same
measurement on the same data. The last three rows have no probe at those shapes
and are `AbsoluteTiming` minimum-of-five.

| probe | before | vs NumPy then | after | vs NumPy now |
|---|---:|---:|---:|---:|
| `pseudoinverse` | **did not finish in 180 s** | — | 19.8 ms | **1/1.07×** |
| `matrixpower` | 14.7 s | 20924× | 827 µs | 1.08× |
| `leastsquares` | **did not finish** | — | 70.9 ms | 1.69× |
| `extract` | 99.7 ms | 97139× | **1.0 µs** | **1/1.35×** |
| `conjugatetranspose` | 430 ms | 195× | 2.18 ms | **1/1.15×** |
| `subdivide` | 1.87 s | 1740× | 707 µs | **1/1.29×** |
| `diagonalmatrix` | 70.2 ms | 320× | 263 µs | **1.09×** |
| `identitymatrix` | 77.3 ms | 334× | 223 µs | 1.16× |
| `unitvector` | 130 ms | 487× | 720 µs | 3.33× |
| `Rescale[v]`, 10⁶ | 2.17 s | — | 3.1 ms | — |
| `ProductLog[v]`, 10⁵ | 2.77 s | — | 14.3 ms | — |

**Every probed row is now within 1.2× of NumPy or ahead of it**, except
`unitvector` at 3.33× on a 720 µs absolute, and `Extract` moved by four orders
of magnitude.

### What the plan got wrong, which is the useful part

- **"They are not packed-aware."** All four of `Extract`, `MatrixPower`,
  `PseudoInverse` and `LeastSquares` were *already* on `AWARE`. They declined
  the buffer themselves — three by `linalg_delist_and_reeval`, one by returning
  `NULL` and letting the **post**-gate materialise. `check_packed_aware.py`
  reports clean on all four; only `MATHILDA_PACK_DIAG=gate` sees it.
- **"C.11: something is not terminating."** Nothing hangs. `PseudoInverse` and
  `LeastSquares` run the *exact* pipeline on 90000 doubles, which is a wrong
  algorithm rather than a loop. `MatrixExp` does not exist and returns in 2 µs;
  the 180 s was the harness's Python side.
- **"C.3: probably an MPFR fallback per element."** `ProductLog` was a **wrong
  answer** — `ProductLog[1.01]` gave −338.392 — and the slow array was the
  *symptom* of the kernel failing and abandoning the buffer. The four Bessel
  heads are blocked by a third thing again: DownValues, which switch off the
  gate's aware test entirely.
- **"C.5: mechanical, no open question."** `Rescale` was not a buffer problem;
  `DiagonalMatrix` of a `Real` diagonal *cannot* be packed.

### Deliberately not done

- **`Array` / `Table` with an exact iterator.** The cost is one interpreter
  evaluation per element, so a buffer alone does not fix it — the body has to
  compile. Compiling with a `CT_REAL` iterator would answer `Table[i^2, …]` in
  floats. The sound version declares the iterator `CT_INT` and accepts only a
  `CT_INT` result type; `CT_REAL` out of `CT_INT` in must be refused (`i/2` is a
  Rational, `Sqrt[i]` is symbolic). Scoped in the HPC plan, not attempted.
- **The Bessel heads.** Making them aware without also giving the *scalar* path
  the libm kernel would leave the packed and plain paths disagreeing in the last
  ulp. The fix is real and needs its own accuracy comparison against MPFR.
- ~~**`DiagonalMatrix` of a `Real` diagonal.**~~ **Done, and it was a wrong
  answer.** I claimed its exact zeros were Mathematica's answer too and filed
  the row under plan item 10.1. Mathematica gives all `Real`s. The rule is
  machine-real contagion — one `Real` makes the whole result `Real`, invented
  zeros included — and `Subdivide` had the identical defect at its endpoints.
  Correcting the exactness made both one dtype, so both now pack: 320× → 1.09×.
  See `tasks/lessons.md`.

### Verification

1. **Differential sweep** — `test_seventh_round_fast_paths` in
   `tests/test_packed_list.c`: 60+ expressions evaluated packed and unpacked
   with byte-identical printed output, including every case that must *decline*
   (mixed-head `Subdivide`, `Real` `DiagonalMatrix`, out-of-range `Extract`,
   symbolic `MatrixPower`, a bad `ConjugateTranspose` spec).
2. **`Subdivide` against `numpy.linspace`** — bit-identical on four intervals
   including a descending one, scaled to exact Integers to see the last bit.
3. **`ProductLog`** — zero residual on `w e^w = x` across ~1100 probe points
   from the branch point at −1/e out to 10³⁰⁰.
4. **Full suite, 396 binaries, rebuilt against the final tree** — one stable
   failure, the documented-stale `simplify_tests` printing expectation
   (`x^2^(3/2)` vs `(x^2)^(3/2)`).
   `moebiusmu_tests` also fails intermittently. **Measured at 4 failures in 12
   runs on an idle machine** — my first note said "load-dependent, 3/3 idle",
   which was three samples and wrong. The cause is `MoebiusMu[10^50 + 1]`
   needing `FactorInteger` to crack a 29-digit composite, and ECM's curve
   search is randomized per process against a fixed budget: the same call
   succeeds 4/4 from a standalone `./Mathilda`. Pre-existing and unrelated —
   no number-theory, factoring or RNG file is touched by this work. Worth its
   own fix (a deterministic seed for the test, or a larger budget), which is
   not this task.
5. **Gates** — `bench_pack`, `bench_eval`, `bench_assoc`,
   `bench_ndarray_linalg`, `make check-c99`, `make check-packed-aware` all pass.
   `ConjugateTranspose` came *out* of the audit's `EXEMPT` table.
6. **Valgrind** — definitely-lost byte-identical to the macOS start-up baseline
   (13,376 B / 418 blocks) over a script exercising every new path.

---

# Follow-up: no two-headed result from a fast path (2026-08-01)

**The invariant**, as stated by the user: a routine handed a packed array /
NDArray returns a scalar or an NDArray — never a mixture of exact and inexact
element heads. A mixture is wrong twice: it disagrees with Mathematica's numeric
tower, and because no uniform buffer holds it, every consumer downstream falls
off the fast path.

`DiagonalMatrix` and `Subdivide` were two instances, fixed above. This is the
sweep for the rest.

**Why it needs a new tool.** The defect is invisible to all three existing
checks. `check_packed_aware.py` asks whether a head opted in — these had. The
differential tests in `test_packed_list.c` compare the packed and unpacked paths
*against each other*, and a producer that invents an exact zero does so on both,
so they agree on the wrong answer. No timing catches it either. What is left is
to look at the element heads, which is what `tools/check_array_exactness.py`
does.

- [x] Write `tools/check_array_exactness.py` — reuses `numeric_sweep`'s 283
      probes for valid call sites, plus a 60-probe corpus chosen to be *read*
      rather than timed (constructors, every matrix-rebuilding linalg routine,
      fills and pads, the narrowing kernels)
- [x] Run it; triage every MIXED row against `wolframscript`
- [x] Fix the divergences; EXEMPT the genuine ones **with the Mathematica output
      as evidence** — the whole reason this file exists is that "Mathematica
      does it too" was once asserted without being checked
- [x] Regression test (`test_no_two_headed_results`) + `make
      check-array-exactness` + documented in `performance.md` §14, the
      changelog, the plan's 10.6/10.6a and `tasks/lessons.md`

## Result

| of 342 probes | before | after |
|---|---:|---:|
| packed | 180 | **184** |
| single-headed plain List | 147 | **148** |
| **MIXED** | **7** | **0** |
| mixed and exempt, with evidence | 8 | 10 |
| no result | 1 | **0** |

(The `no result` row was `NullSpace[A300]`, which turned out to be the same
`pack_offer` breakage rather than a timeout.)

Six routines were inventing an exact element inside a machine-real result:

| | was | now | Mathematica |
|---|---|---|---|
| `VandermondeMatrix[{1., 2., 3.}]` | `x^0` column exact | all `Real` | agrees |
| `HankelMatrix[{1, 2, 3.}]` | pad zeros exact | all `Real` | agrees |
| `ToeplitzMatrix[{1, 2, 3.}]` | entries not widened | all `Real` | agrees |
| `MatrixPower[{{2., 0.}, {0., 2.}}, 0]` | `{{1, 0}, {0, 1}}` | all `Real` | agrees |
| `NullSpace[{{1., 2.}, {2., 4.}}]` | `{{-2., 1}}` | `{{-2., 1.}}` | agrees |
| `RowReduce[{{2., 0.}, {0., 4.}}]` | `{{1, 0.}, {0., 1}}` | all `Real` | **differs** |

All six now pack. The rule — invented elements take the input's exactness — is
`common_has_machine_real` / `common_machine_real_value` in `src/common.h`;
`inverse_onestep` had been doing it by hand since the fifth sweep.

**`RowReduce` is a stated divergence.** Mathematica writes the exact `1` and `0`
regardless: its `RowReduce[{{2., 4.}, {1., 3.}}]` has heads
`{Integer, Real, Integer, Integer}` and `PackedArrayQ` `False`, so its own RREF
of a machine matrix is two-headed and unpacked. Mathilda follows the project
rule instead. Values agree; only the structural entries' heads differ.

**The eight exemptions are where the rule stops**, each carrying the Mathematica
output in the tool: `Chop`'s zero is exact on purpose; the
`PadRight`/`PadLeft`/`ArrayPad`/`Riffle`/`Insert`/`Append`/`ReplacePart` family
places a caller-supplied element verbatim; `Join` mixes what it was given;
`Tally` returns `{value, count}` pairs; `LUDecomposition`'s exact elements are
its permutation vector.

**Still open:** `NullSpace` of an inexact matrix is not orthonormalized
(Mathematica gives `{{-0.894…, 0.447…}}` where we give `{{-2., 1.}}`) — same
basis, different normalisation, a larger change than exactness. And
`CholeskyDecomposition`, `Orthogonalize` and `ArrayPad` return **unevaluated**
on a machine matrix; that is the C.12 coverage gap, not an exactness one.
