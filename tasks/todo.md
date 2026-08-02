# Ninth round — the 26 heads that declined a visible NDArray (2026-08-02)

The eighth round's record is below, and `docs/design/performance.md` §15; earlier
rounds are §8–§14 and `docs/experiments/`.

## Why this round exists

The eighth round's three-surface audit ended with a register of six finding
kinds. The largest correctness item left open was `ND-UNSUPPORTED`: **26 heads
that answer on a plain List and on a packed List, and leave the call
unevaluated on a visible `NDArray`.**

The eighth round wrote down a general fix and declined to take it:

> `eval.c` already has a post-gate that materialises a *packed* argument when a
> node comes to rest… The same post-gate for a *visible* array would close all
> 26 at once. The hazard is that `NDArray[...]` is itself a head that comes to
> rest holding an array… That is a change to the core evaluator and deserves
> its own round.

**That is not the fix taken here either, and the reason is not the hazard.** A
post-gate makes the call *evaluate*; it does not make it *fast*. It would have
turned 26 unevaluated calls into 26 calls that materialise 10⁶ `Expr` nodes and
run the generic List code — closing the audit row and leaving the cost exactly
where it was. Every head below got a real buffer path instead.

## Plan

- [x] Batch A — the engine gap: binary kernels accepted only ONE array operand
- [x] Batch B — the integer domain: GCD, LCM, DivisorSigma, EulerPhi, MoebiusMu,
      IntegerLength, PowerMod, Prime, IntegerDigits
- [x] Batch C — the sign predicates: Positive, Negative, NonNegative, NonPositive
- [x] Batch D — structural: Ratios, Append, Prepend, Catenate, TakeLargest,
      TakeSmallest, Counts, Inner, RandomSample, RandomChoice
- [x] Batch E — hypergeometric: `Hypergeometric{0F1,1F1,2F1}` via `PFQ`
- [x] Both surfaces: `AWARE` / `INT64_OK` in `pack.c`, so a packed List reaches
      the same path a visible array does
- [x] Differential test of every head on all three surfaces
- [x] Regression tests, docs + changelog

---

## What was built

### Batch A — one array operand was hard-coded into the engine

`ndarray_map_binary` required **exactly one** array and one broadcast scalar,
and `eval.c`'s dispatch enforced the same with an XOR. So a kernel registered
for a genuinely two-argument function was unreachable whenever both arguments
were arrays — which is the ordinary way to call one. `ArcTan[v, w]` is numpy's
`arctan2` of two vectors and `Beta[p, q]` takes two; both fell straight through
to the unevaluated call, with `NDKB_ArcTan` and `NDKB_Beta` registered and
unused. **Fifteen registered kernels were in that position**, so this was one
engine gap rather than fifteen head-level ones.

`ndarray_map_binary2` mirrors the scalar-broadcast form branch for branch —
exact-integer arm, int64 guard, complex, real-closed, escaping-real — reading
the second operand at the same index instead of hoisting it out of the loop,
and promoting the two input dtypes for the result.

A second, smaller gap fell out of it: the `!k->cplx` test at the top meant
"sentinel, degrade", which made an **integer-only** kernel inexpressible. It
moved down to the branches that actually call `cplx`, unchanged for every
kernel that has one, and Batch B's six kernels became possible.

### Batch B — the integer domain had no kernel at all

New translation unit `src/ndinteger.c`. `ndkernels.c` owns the
libc-expressible `double` set and says so at the top; these are the opposite
kind — defined on ℤ, answers are exact Integers, and a `double` intermediate is
not an approximation of the answer but a different answer. Written in int64 end
to end through `ci_*_i64`, abandoning the whole array on overflow so the List
path and GMP still answer exactly.

Factorisation is trial division rather than GMP: the right choice for one
number of unbounded size is the wrong one for 200 000 of them, where the
per-call `mpz_init`/`clear` alone dominates. Capped at 10¹² (10⁶ divisions);
past it the kernel declines and GMP answers, so the ceiling costs speed and
never correctness.

Three heads cannot be an element kernel and got their own entry point:

| | why not a kernel |
|---|---|
| `Prime` | per-element `Prime[n]` is a Meissel count each time. Over an array the answer is a **sieve** to Rosser's bound for the largest index, then a gather — a different algorithm, not a different loop |
| `PowerMod` | ternary; `NDUnary`/`NDBinaryKernel` cannot express three arguments |
| `IntegerDigits` | the result is **ragged**, so no buffer holds it. The win is on the input side |

### Batch C — the sign predicates

`Positive`/`Negative`/`NonNegative`/`NonPositive` answer with a List of
True/False, which no dtype holds (there is no boolean dtype — §13 gap C.1). The
fast path is therefore one-sided: read the buffer directly, no `Expr` per input,
no evaluator round-trip. Two shared symbol nodes carry the whole output by
reference, since `expr_copy` is a refcount bump — 10⁶ increments instead of 10⁶
allocations.

### Batch D — structural

`Ratios` divides in place (integer dtypes decline: the ratios are exact
Rationals). `Append`/`Prepend` are one allocation and two memcpys — they were
named in `pack.c`'s "correct by omission" list, which is right for correctness
and costs a full 10⁶-element boxing to add one element. `Catenate` is a reshape
for a rank-2 buffer and a `Join` for a List of arrays, which is the same
operation seen from the two surfaces. `TakeLargest`/`TakeSmallest` keep a
bounded heap, O(n log k) — ahead of the full sort NumPy's own idiom for this row
uses. `Counts` relabels `Tally`'s machine-word count. `Inner[Times, a, b, Plus]`
**is** a Dot, so it is one.

`RandomSample`/`RandomChoice` gather from the buffer **through the same draw
sequence the List path uses** — `fisher_yates_sample`, `random_index` — not a
new shuffle. A producer that consumed the generator differently would make
`SeedRandom[1]; RandomSample[v]` give one permutation for a packed `v` and
another for the identical plain `v`: a surface disagreement in the one place a
user cannot see it coming. Verified seeded, on all three surfaces.

### Batch E — hypergeometric

`Hypergeometric{0F1,1F1,2F1}` all rewrite to `HypergeometricPFQ` before anything
else sees them, so the path belongs in `PFQ`. **Every parameter must be
inexact**, and that is the correctness argument rather than a convenience:
`try_cancel`, `try_terminate` and `try_reduce` are not slower routes to the same
number — they answer with exact or closed forms (`2F1[1,1,2,z]` is
`-Log[1-z]/z`) and fire on exact parameters. A Real parameter cannot match them,
so the series really is what the scalar path would have computed. Checked both
ways, including the terminating and closed-form cases.

---

## Three defects the differential test found, which the audit would not have

Each was introduced by this round's own work and each is the same shape: giving
a head a fast path changes what reaches its *slow* path too.

1. **`Counts` stayed unevaluated on exactly the data it exists for.** `Tally`
   answers in two shapes — over an int64 buffer its `{key, count}` pairs are
   themselves machine words, so it hands back a rank-2 **array**, while the
   float64 route builds a List of pairs. Reading only the List shape left the
   integer case dead.
2. **`Append[NDArray[{1., 2.}], 0]` coerced the exact `0` to `0.`.** The decline
   path used `ndstruct_delist_repack`, which re-packs a visible array's result
   at the source dtype — and the whole reason it declined is that the appended
   element does not belong to that dtype. Mixed exact/inexact is the correct
   answer and no buffer holds it: `ndarray_delist_and_reeval`, not the repack.
3. **`Inner` with any other operator pair fell out unevaluated.** Marking a head
   `AWARE` stops the gate materialising, so the buffer now reaches code that
   tests `type == EXPR_FUNCTION`. **This is the standing hazard of the whole
   round** and it applies to every head added to `AWARE`: `Prime`, `PowerMod`,
   `IntegerDigits` and `HypergeometricPFQ` each needed the same explicit
   degrade, written before the test found them because the pattern was already
   known from `Inner`.

### And a fourth head, found by checking whether the list was right

The eighth round listed the 6×6 `Dot`/`Inverse`/`LinearSolve` rows among the 26.
Re-probed directly, all three answer a visible `NDArray` correctly and always
did — so the row was wrong, and running down *why* found a real bug in a head
nobody had suspected. The probe's checksum is `N[Total[Flatten[{r}]]]` over a
`List` of 200 result arrays, and **`Flatten` treated a visible `NDArray` nested
in an ordinary `List` as an atom**:

```
Flatten[{{1., 2.}, {3., 4.}}]                    -> {1., 2., 3., 4.}
Flatten[{NDArray[{1., 2.}], NDArray[{3., 4.}]}]  -> unchanged
```

`head_is` is false for an `EXPR_NDARRAY`. It is a list of values by every other
measure (`ArrayQ`, `Dimensions`, `Length` all say so), so `flatten_rec` descends
into it now. The packed form was never affected — the no-nesting invariant keeps
a packed node out of a plain `List` — which is **the same asymmetry as the eighth
round's int64 kernels and pattern family, in a third subsystem**. That is three
independent instances of one shape, and the argument for keeping the audit
pointed at surfaces rather than sites.

## Measured

`tools/nd_surface_audit.py`, all three representations, Apple M-series — the
tool's own probes and its own gain/skew arithmetic, not a hand-rolled harness
(see the correction note at the end of this file for why that distinction cost
me two wrong rows).

| probe | plain | packed | visible | gain | skew |
|---|---:|---:|---:|---:|---:|
| `inner_dot` | 243 µs | **1 µs** | 2 µs | 243× | |
| `beta_fn` | 1.299 s | **4.87 ms** | 5.10 ms | 267× | |
| `append` | 96.99 ms | **702 µs** | 732 µs | 138× | |
| `ratios` | 272.2 ms | 2.08 ms | **1.62 ms** | 131× | |
| `arctan2` | 311.2 ms | **2.59 ms** | 3.39 ms | 120× | |
| `integerlength` | 92.68 ms | **855 µs** | 1.01 ms | 108× | |
| `prepend` | 109.2 ms | **1.44 ms** | 2.00 ms | 76× | |
| `catenate` | 233.6 ms | 5.62 ms | **3.43 ms** | 42× | 1.64× `SKEW` |
| `divisorsigma` | 118.5 ms | **3.30 ms** | 3.30 ms | 36× | |
| `counts` | 41.05 ms | 1.88 ms | **1.83 ms** | 22× | |
| `randomchoice` | 23.51 ms | **1.24 ms** | 1.25 ms | 19× | |
| `takesmallest` | 9.82 ms | **556 µs** | 634 µs | 18× | |
| `randomsample` | 26.04 ms | 1.64 ms | **1.60 ms** | 16× | |
| `powermod` | 109.9 ms | 7.15 ms | **6.94 ms** | 15× | |
| `hyper1f1` | 856.0 ms | 67.56 ms | **66.88 ms** | 12.7× | |
| `hyper2f1` | 1.528 s | **120.6 ms** | 124.3 ms | 12.7× | |
| `prime_arr` | 7.52 ms | 620 µs | **580 µs** | 12.1× | |
| `moebiusmu` | 35.56 ms | **3.09 ms** | 3.47 ms | 11.5× | |
| `eulerphi` | 26.07 ms | **3.64 ms** | 4.01 ms | 7.2× | |
| `takelargest` | 10.60 ms | **2.50 ms** | 2.66 ms | 4.2× | |
| `positive` | 275.1 ms | **66.42 ms** | 72.18 ms | 4.1× | |
| `nonnegative` | 272.9 ms | **66.68 ms** | 73.82 ms | 4.1× | |
| `negative` | 271.6 ms | **67.28 ms** | 73.14 ms | 4.0× | |
| `lcm_arr` | 72.53 ms | 29.35 ms | **8.30 ms** | 2.5× | 3.54× `SKEW` |
| `integerdigits` | 26.53 ms | **10.81 ms** | 11.41 ms | 2.5× | |
| `gcd_arr` | 81.42 ms | 38.60 ms | **16.74 ms** | 2.1× | 2.31× `SKEW` |

`gain` is the better machine surface against the plain-List floor, as the tool
computes it. Only **four** rows sit under 5×, and both reasons are structural
rather than unfinished work:

* `positive`/`negative`/`nonnegative` (≈4×) and `integerdigits` (2.5×) must
  build one `Expr` per **output** element whatever happens — a list of
  `True`/`False` and a ragged list of digit lists are not buffers. Reading the
  input off the buffer is the whole of the available win.
* `gcd_arr`/`lcm_arr` are held down by their *packed* column only; on the
  visible surface they are 4.9× and 8.7×. That gap is the `Orderless` finding
  recorded below.

A float64 hot lane was added to both binary map chunks (the twin of
`ndu_hot_chunk`). A/B'd by building the tree with and without it: `ArcTan[v, w]`
3.26 ms → 2.49 ms (23%), `Beta[p, p]` 5.27 ms → 4.94 ms (6%). Beta gains less
because three `lgamma`s per element dominate what the marshalling costs, where
`atan2` does not — so the lane is worth most exactly where the kernel is
cheapest.

### A 1-ulp divergence found and closed on the way

`ndk_ArcTan2_c` computed `arg` via `csqrt`/`clog` where the scalar builtin calls
`atan2`, and landed 1–2 ulp away on real data — **68 of 400 elements** differed
on the sweep's `arctan2` probe. Under the sweep's 1e-5 tolerance it never
showed, and it predates this round: the array/scalar route had it too (44 of
400). The kernel now calls `atan2` for real inputs, which is what this file's
own contract says a kernel does, and both routes are bit-identical to the scalar
path. `atan2` also gets `(0, 0)`, the negative-zero convention and the
infinities right by definition, where the old form declined on `s == 0` and
abandoned the whole array for one such pair.

`Beta` and `1F1`/`0F1` remain ~1 ulp from their MPFR twins (2.1e-16 and 1.6e-16
relative). That is the register's "MPFR-vs-libm accuracy comparison", now done:
it is the same relationship `Gamma`, `Erf` and the Bessel kernels have always
had, and two orders inside the sweep's tolerance.

## Both open skews, closed

### 1. `Orderless` materialised the buffer to compare it against a scalar

`GCD` 2.3× and `LCM` 3.5× slower packed than visible. `Mod`/`Quotient` — same
kernels, same shapes, **not** `Orderless` — showed no gap, and that control is
what named the cause in one step.

`expr_compare` (`src/sort.c` step 0) orders a packed list as the List it is,
correctly, but reached that via `ndarray_to_nested_list` when only one operand
was packed. The evaluator sorts `GCD[1234, cv]` before dispatch, so 200 000
`Expr` nodes were built and discarded per call — to settle an ordering that
steps 1–2 decide by tier membership, without reading a single element.

Fixed by comparing against an **empty `List` stand-in**, and only where the
answer provably cannot depend on the elements (other operand a number, a literal
`Complex`, or a `String`). A hardcoded `return 1` would have been shorter and
would have restated an ordering rule in a second place; the stand-in lands in
the same tier and lets the existing code compute it. A `Symbol` or general
expression on the other side reaches step 3's polynomial-degree walk, which does
read elements, and still materialises.

| n = 200 000 | before | after | visible |
|---|---:|---:|---:|
| `GCD[cv, 1234]` | 46.7 ms | **17.1 ms** | 18.6 ms |
| `LCM[cv, 12]` | 31.5 ms | **7.44 ms** | 7.48 ms |

Ordering verified unchanged against the materialising path across every tier —
`Integer`, negative, `Rational`, `Real`, exact and inexact `Complex`, `String`,
`Symbol`, a general expression, a plain `List`, a second packed list — both
directions, plus `Sort` of a mixed bag and `Orderless` through `Plus`. Identical
with packing on and off.

### 2. `PACK_MIN_ELEMENTS` 250 → 4

A 6×6 is 36 elements, under the threshold, so it never packed and LAPACK was
unreachable from it. The old value's comment gave its own reasoning — chosen for
*blast radius*, not cost, break-even already known to be around n = 2 — and that
is what proved wrong: the margin was deliberate but measured in neither
direction.

| `Table[…, {200}]` over a 6×6 | 250 | **4** | speed-up |
|---|---:|---:|---:|
| `Det[A6]` | 102.8 ms | **0.189 ms** | 544× |
| `Inverse[A6]` | 120.8 ms | **0.473 ms** | 255× |
| `LinearSolve[A6, b6]` | 134.2 ms | **0.397 ms** | 338× |
| `A6 . A6` | 12.9 ms | **0.108 ms** | 120× |

The benefit side was then swept rather than assumed: 250, 64, 36, 16, 8, 4, 2
against pattern matching, rule application, `Cases`/`MatchQ`, `Table`, `Expand`,
`Solve`, `D`, `Integrate`, `Sort`, `Join`, `Nest`, `Counts`, `Simplify`. **Down
to 4, no regression on any**, and the nested `Table` gets 19% faster (64
elements now packs). At 2 the linalg stops improving and `Integrate` gives back
7%. So 4 — the element count of a 2×2 matrix, i.e. *any matrix packs*.

Residual skew under 5% (`Det` 158 µs packed against 166 µs visible), from
118–587×. What little remains is not the threshold:
`Table[…, {200}]` over packed results absorbs the 200 outputs into a rank-3
buffer, which the visible column does not do — useful work, not waste.

`MATHILDA_PACK_MIN` overrides it per session (same family as `MATHILDA_NO_PACK`
/ `MATHILDA_PACK_DIAG`), because this is the one packing parameter whose right
value is empirical and re-sweeping it should not need a rebuild.

### A third change the second one forced: signed zero in machine linalg results

Lowering the threshold put a 2×2 on the LAPACK path, and LAPACK produces `-0.0`
wherever a zero entry came out of a subtraction or a sign flip:

```
In[1]:= Inverse[Table[If[i == j, 2., 0.], {i, 2}, {j, 2}]]
Out[1]= {{0.5, -0.0}, {0.0, 0.5}}       (* packed  *)
Out[1]= {{0.5,  0.0}, {0.0, 0.5}}       (* plain   *)
```

The sign of a zero in a matrix result carries no information — the entry is
exactly zero and which way it was reached is an artefact of the pivoting — but
it **prints**, and it printed differently from the exact path answering the same
expression. That is precisely what `pack.h`'s "a representation may never change
a value" exists to prevent.

The behaviour is **pre-existing**, not introduced here: a 20×20 diagonal matrix
gave `-0.0` at the old threshold too, because it was already over it. What
changed is that the sizes where it shows are now the sizes people read. So
`na_build_vector` and `na_build_matrix` — the two shared builders behind every
machine linalg result — normalise it, which fixes large matrices as well.

Found by a regression test written for the threshold change, not by the
benchmark: the value is unchanged under `SameQ`, so nothing that compares
results would have caught it.

### Still open: a literal list never packs, at any size

A **different** gate, found while checking the above. `NDArrayQ[Table[1., {400}]]`
is `True`; `NDArrayQ[{1., 2., 3., 4., 5., 6.}]` — or a 400-element literal — is
`False`. Packing is opt-in per **producer**; a `List` from the parser has none.
So `Inverse[{{4., 1.}, {1., 3.}}]` typed literally still takes the unpacked path
however low the threshold goes, and this change helps only *computed* small
matrices. Closing it means offering every list node that comes to rest, which is
the blast radius the producer-opt-in design exists to avoid.

### The audit's verdict

`tools/nd_surface_audit.py --only …` over the 30 probes for the heads changed
here plus four controls, all three representations, after every fix in this
round:

| kind | before | after |
|---|---:|---:|
| `ND-UNSUPPORTED` | 26 | **0** |
| `SKEW` | 7 | **0** |
| `NO-PATH` | 4 | **0** |
| `DISAGREE` | 0 | **0** |
| `NO-ANSWER` | 0 | **0** |
| `ND-SLOW` | 0 | 4 — noise, see below |

The small-matrix rows are the clearest single result: `Det[A6]` **617×**,
`LinearSolve[A6, b6]` **341×**, `Inverse[A6]` **277×**, `A6 . A6` **119×**
against the plain-List floor, with packed and visible now within 5% of each
other (158 µs against 166 µs for `Det`). `GCD` and `LCM` come in at 16.4/16.6 ms
and 7.46/7.44 ms — the two surfaces level, where they were 2.3× and 3.5× apart.

**The four `ND-SLOW` rows are a measurement artefact and were checked rather
than reported.** `nonnegative` tripped at 0.59× while `positive` and `negative`
— the same operation, the same data — sat at 0.74 and 0.73, which is not how a
real property behaves. That run shared the machine with a
`check-array-exactness` pass: re-measured idle, `nonnegative` costs 24.6 ms
rather than the 65–110 ms the audit recorded, and all four probes come out level
across two rounds:

| 10⁶ elements, idle | packed | visible |
|---|---:|---:|
| `Ratios[v]` | 1.65 / 2.08 ms | 1.93 / 1.98 ms |
| `Counts[jv]` | 1.86 / 2.02 ms | 1.85 / 2.01 ms |
| `NonNegative[u]` | 24.6 / 28.5 ms | 24.8 / 26.7 ms |
| `Union[kv]` | 4.54 / 4.27 ms | 4.12 / 3.84 ms |

A timing tool sharing a machine with anything else is measuring the other thing
too — the same lesson as the hand-rolled harness earlier in this round, one
level up.

The full 284-probe three-surface run does not finish inside 90 minutes (killed
at 228/284), so the targeted subset is what is recorded. Heads outside this
round's set rest on the eighth round's run.

## Verification

* Every head above differential-tested against the List surface on **all three
  representations** — plain, packed, visible — including the negative, zero,
  overflow, bignum-escape, wrong-dtype, ragged and seeded-RNG cases.
* `make check-c99` — pass.
* `make check-packed-aware` — pass; `AWARE=157`, `INT64_OK=96`.
* `make check-array-exactness` — 342 probes, **0 MIXED**.
* `make check-nd-surfaces` — see the run recorded in
  `docs/design/performance.md` §16.
* Full test suite — see below.


## Correction — the first measurement pass was wrong, and how

The `ArcTan[v, w]` and `Beta[p, p]` rows were first reported at **3.6× and
2.7×**. They are **112× and 222×**. The harness built the second `NDArray`
operand *inside* the timed expression:

```
t["ArcTan visible", ArcTan[vn, NDArray[Reverse[v], DataType -> "float64"]]]
```

`Reverse` over a plain 10⁶ List plus the pack is ~90 ms; the kernel is ~2.5 ms.
The measurement was of the conversion, and every ratio taken from it was a
statement about `Reverse`.

The wrong numbers were the smaller problem. They supported a confident and
wrong **conclusion** — "`Beta` and `ArcTan` are libm-bound, the marshalling is
not the cost" — which was then written into `ndarray.c` as a comment, into the
changelog as a finding, and into `performance.md` as guidance for the next
round. It also made a float64 hot lane in the binary map chunks look worthless:
measured 8% and nothing, where a real A/B (building the tree with and without
the lane) gives `ArcTan` 3.26 → 2.49 ms (23%) and `Beta` 5.27 → 4.94 ms (6%).

Caught by `nd_surface_audit.py` disagreeing by two orders of magnitude on the
same probes. It builds every array in a preamble, which is precisely the thing
the hand-rolled harness got wrong — so the tool's table replaced mine
everywhere, and the rule is: **hoist every conversion, and when the project's
tool disagrees with your harness, the tool is right until proven otherwise.**

---
---

# Eighth round — the packed / NDArray surface, audited rather than assumed (2026-08-01)

The seventh round's record is `docs/design/performance.md` §14; earlier rounds
are §8–§13 and `docs/experiments/`.

## Why this round exists

`DeleteDuplicates` sat at 72× NumPy through four sweeps with **no NDArray path
at all**. Neither existing audit could see it:

* `tools/check_packed_aware.py` reads *dispatch sites* out of the source, so a
  head with no fast path on either surface has nothing to read and is invisible
  to it by construction.
* `tools/numeric_coverage.py` joins the *registries*. Registration is not speed.

Both are static. The question they cannot answer is the one that matters: **run
the same expression on each representation and see whether it is actually fast.**

## The three representations, and why they can disagree

| | gate behaviour | consequence |
|---|---|---|
| plain `List` | — | one `Expr` per element, the floor |
| **packed** `List` | `src/eval.c` **materialises** it for any head not on `pack.c`'s `AWARE` list | a missing opt-in fails *safe*, therefore *silent* |
| **visible `NDArray`** | the gate never touches it (`is_packed_list` is false) | reaches the builtin whatever `AWARE` says |

Opposite signs — so a head can be fast on one surface and slow on the other on
identical data. Not hypothetical: `Fourier` had a complete NDArray path, was
absent from `AWARE`, and ran 15.7× slower packed.

## Plan

- [x] `tools/nd_surface_audit.py` — the same 284 probes as `numeric_sweep.py`,
      on all three surfaces, flagging `NO-PATH`, `SKEW`, `ND-SLOW`, `DISAGREE`
- [x] `--survival`: does the head *return* a packed array, or hand back a plain
      List that `ToNDArray` would have packed? Packing is a chain; a producer
      that drops it makes every consumer slow, not itself
- [x] Run the full three-surface audit; triage every finding
- [x] Fix the top three, in the order the measurements put them
- [x] Regression tests, a `make` target, docs + changelog

## Findings

### 1. The integer pipeline is unpacked from the first `Mod`

`jv = Mod[Range[10^6]*7919, 1000]` comes back **plain**. `Mod` is
`packed_aware` (it has `NDKB_Mod`) but is not in `INT64_OK` — rightly, because
its kernel computes `m - n*floor(m/n)` in `double` and would write `{0., 1.}`
where the List gives the exact `{0, 1}`. So the gate materialises, the result is
10⁶ boxed Integers, and every consumer downstream pays:

| | packed | visible NDArray | skew |
|---|---:|---:|---:|
| `Union[kv]` | 807 ms | 4.44 ms | **182×** |
| `Tally[jv]` | 68.1 ms | 1.49 ms | **45.8×** |
| `DeleteDuplicates[jv]` | 67.1 ms | 1.98 ms | **34.0×** |

The set-op work of 2026-08-01 is not at fault: those heads are on `AWARE` *and*
`INT64_OK` and are fast when handed a buffer. They never are.

### 2. `--survival`, elementwise + integer: 13 producers drop packing

`Im`, `Boole`, `ArcTan[v, w]`, `GCD`, `LCM`, `Mod`, `Quotient`,
`IntegerLength`, `PowerMod`, `EulerPhi`, `MoebiusMu`, `DivisorSigma`, `Prime`.

### 3. Not an attribute bug — 96 numeric heads do not exist

`BitAnd`/`BitOr`/`BitXor`/`BitNot`/`BitShiftLeft`/`BitShiftRight`, `CubeRoot`,
`Log2`, `Log10`, `Surd`, `Unitize`, `LogisticSigmoid`, `Gudermannian`,
`Haversine`, `KroneckerDelta`, `Threshold`, `Quantile`, `Correlation`,
`CholeskyDecomposition`, `MatrixExp`, `Ordering`, … return unevaluated because
the symbol is undefined, which is why `Attributes[…]` reads `{}`. Their sweep
rows are the evaluator declining to answer, not slow code.

This is `performance.md` §13 register item 6 — a **coverage** gap rather than a
packing one. Recorded here, not closed here.

---

## What was fixed

### 1. Every real kernel truncated a visible int64 `NDArray` — a silent wrong answer

```
Sin[NDArray[{1, 2, 3}, DataType -> "int64"]]  ->  NDArray[{0, 0, 0}]
Exp[NDArray[{1, 2, 3}, DataType -> "int64"]]  ->  NDArray[{2, 7, 20}]
```

`ndarray_map_unary`/`_binary` size the output from the **input** dtype and write
through `ndt_set`, whose `NDT_INT64` case is `(int64_t)re`. Every `real_closed`
kernel without an exact integer arm — 56 of the 85 — truncated. `Cos`, `Tanh`,
`Erf`, `BesselJ`,
`Log[b, ·]` all did.

The packed surface was never exposed to it, because the gate materialises an
int64 buffer for any head without `packed_int64_ok`. **The guard that makes the
packed surface safe is exactly why the visible one was unguarded.** Both map
functions now decline an int64 input with no exact arm.

### 2. `Mod` was not slow — it made everything downstream slow

| 10⁶ elements | before | after |
|---|---:|---:|
| `Union[kv]` | 807 ms | **4.58 ms** |
| `Tally[jv]` | 68.1 ms | **1.89 ms** |
| `DeleteDuplicates[jv]` | 67.1 ms | **1.80 ms** |

`NDBinaryKernel` gained the exact-integer arms `NDUnaryKernel` has had since the
narrowing kernels. `Mod` takes `to_int_i` only (its exactness follows the
argument); `Quotient` takes both (it narrows to an exact Integer whatever it is
given) and came **off** `NOT_AWARE`, closing the threshold divergence that list
existed for — `Quotient[Range[1., 300.], 2]` and `Quotient[Range[1., 200.], 2]`
now agree, both `{0, 1, 1, 2}`.

### 3. The set operations were unreachable from a visible `NDArray`

`setop_i64` tested `is_packed_list`. `Union` 850 ms → **8.47 ms**,
`DeleteDuplicates` 147 ms → **2.34 ms**. A *dispatch* gate wants `is_ndarray`;
only a *presentation* decision wants `is_packed_list`.

## The register — 38 producers that still drop packing

From `make check-nd-surfaces` (`--survival`), filtered to results at or over the
250-element threshold. Grouped by cause, since the fix differs by group:

| group | probes | note |
|---|---|---|
| **structural, has a packed sibling** | `catenate`, `ratios`, `split`, `position`, `gatherby`, `sortby`, `mapindexed`, `thread_plus` | `Join` packs and `Catenate` does not; `Differences` packs and `Ratios` does not. **Checked, not assumed:** `{v, w}` is absorbed into a rank-2 packed array by `pack_sniff`, so `Catenate` receives an *NDArray* where its code tests `type == EXPR_FUNCTION` — for a packed argument it is exactly `Flatten[arr, 1]`, a memcpy. `SortBy` is already a documented `EXEMPT` in `check_packed_aware.py` (its NDArray path does not sort) |
| **insert/delete family** | `replacepart`, `insert`, `delete`, `append`, `prepend` | register C.6, 153–198× — "correct by omission" in `pack.c` |
| **number theory over int64** | `gcd_arr`, `lcm_arr`, `integerlength`, `powermod`, `eulerphi`, `moebiusmu`, `divisorsigma`, `prime_arr` | exact integer domain, no kernel at all. The natural next batch after `Mod`/`Quotient` |
| **special functions** | `besselj`, `bessely`, `besseli`, `besselk`, `beta_fn`, `hyper2f1`, `hyper1f1` | register C.3 — blocked by `DownValues` from `internal/*.m`, and needs the MPFR-vs-libm accuracy comparison first |
| **linalg results** | `qr`, `svd`, `eigenvalues`, `eigenvectors` | the decomposition returns `{q, r}` and **neither matrix is packed** — so the outer `List` has nothing for `pack_sniff` to absorb. Not a nesting limitation: `{v, w}` of two packed vectors *does* become a rank-2 array. The routines simply do not pack their outputs |
| **dtype / representation** | `n_of_list` (`N[iv]`, 10⁶), `im`, `boole`, `arctan2` | `N[list]` is a dtype cast; `Im` needs the narrowing treatment `Floor` got; `Boole` needs the boolean dtype (C.1); `ArcTan[v, w]` is two array operands, which `ndarray_map_binary` does not accept |
| **random** | `randomsample`, `randomchoice` | producers without a `pack_offer` |

Six more (`minmax`, `extract`, `takelargest`, `takesmallest`, `cross`, `det6`)
are packable but **under** the threshold, where not packing is correct — the
first run reported them because `ToNDArray` deliberately ignores the threshold.
The tool filters them now.

---

## The audit's own register — 284 probes, three surfaces, after the fixes

Findings re-classified from the completed run's JSON. Four kinds, and a fifth
that the first classifier got wrong.

| kind | n | what it means |
|---|---:|---|
| `DISAGREE` | 11 → **4** | the surfaces compute different things. Seven were the pattern family, fixed below |
| `ND-UNSUPPORTED` | 26 | answers on plain and packed, **unevaluated** on a visible `NDArray` |
| `ND-SLOW` | 9 | visible array materially slower than packed |
| `SKEW` | 54 | packed materially slower than visible — mostly heads with no path on either, where both columns are noise |
| `NO-PATH` | 98 | packed no faster than plain |
| `NO-ANSWER` | 29 | answered on **no** surface: the 96 undefined heads (§13 register C.12) |

### 4. The pattern family returned confident wrong answers on a visible `NDArray`

The severe one, found by `DISAGREE`. The matcher walks `data.function.args`;
an `EXPR_NDARRAY` has none, so each head searched an expression with no elements
and reported success at finding nothing:

| | visible `NDArray` | `List` |
|---|---|---|
| `MemberQ[·, 5.]` | **False** | True |
| `Count[·, 5]` | **0** | 1 |
| `Position[·, 5]` | **{}** | `{{5}}` |
| `Cases[·, 5]` | **{}** | `{5}` |
| `FreeQ[·, 5]` | **True** | False |

Worse than declining — a wrong answer that looks like a right one. The packed
form was never affected: these heads are not on `AWARE`, so the gate
materialises their arguments first. **The same asymmetry as the int64 kernels,
in a completely different subsystem**, which is the argument for auditing the
surfaces rather than the sites. Each now materialises a visible array itself
(`patterns_delist_visible`); `FreeQ` in `funcprog.c` likewise.

### `ND-UNSUPPORTED`: 26 heads decline a visible `NDArray` outright

`Positive`/`Negative`/`NonNegative`, `Counts`, `RandomSample`, `RandomChoice`,
`IntegerDigits`, `TakeLargest`/`TakeSmallest`, `Append`, `Prepend`, `Ratios`,
`Inner`, `Catenate`, `ArcTan[a, b]`, `GCD`, `LCM`, `IntegerLength`, `PowerMod`,
`EulerPhi`, `MoebiusMu`, `DivisorSigma`, `Prime`, `Beta`,
`Hypergeometric1F1`/`2F1`, and the 6×6 `Dot`/`Inverse`/`LinearSolve` rows.

These fail **loudly** (the call stays unevaluated), which is why they rank below
the wrong-answer class.

**The general fix, and why it was not taken here.** `eval.c` already has a
post-gate that materialises a *packed* argument when a node comes to rest, on
the reasoning that the head has now declined the buffer. The same post-gate for
a *visible* array would close all 26 at once. The hazard is that `NDArray[...]`
is itself a head that comes to rest holding an array — its own constructor — so
the rule needs an exclusion list, and getting that wrong breaks the
representation everywhere. That is a change to the core evaluator and deserves
its own round with its own sweep, not the tail of this one.

*(Closed by the ninth round above — with buffer paths rather than a post-gate,
for the reason recorded there.)*

### `ND-SLOW`: 9 heads slower on the visible surface

`max2`, `min2`, `take`, `part_gather`, `diagonalmatrix`, `map_pure`,
`map_lambda`, `mapthread_plus`, `outer_sub`. `Map` is the one with a named
cause: `numloop.c`'s fast path tests `is_packed_list` and so declines a visible
array — the same predicate confusion that cost the set operations 145×.

## Verification

Re-measured after the fixes, all three surfaces, same probes:

| probe | plain | packed | visible `NDArray` | gain | skew | findings |
|---|---:|---:|---:|---:|---:|---|
| `sin` | 372.8 ms | 1.79 ms | 1.66 ms | 208× | 1.08× | — |
| `total` | 85.6 ms | 247 µs | 277 µs | 347× | 0.89× | — |
| `union` | 596.4 ms | 4.90 ms | 5.53 ms | 122× | 0.89× | — |
| `tally` | 54.7 ms | 1.69 ms | 1.93 ms | 32× | 0.87× | — |
| `deleteduplicates` | 51.5 ms | 1.96 ms | 1.92 ms | 26× | 1.02× | — |
| `mod_int` | 597.2 ms | 11.5 ms | 11.7 ms | 52× | 0.98× | — |
| `quotient_int` | 650.1 ms | 12.9 ms | 12.5 ms | 50× | 1.03× | — |
| `memberq` | 112.1 ms | 165.9 ms | 206.5 ms | 0.7× | 0.80× | `NO-PATH` |
| `position` | 193.2 ms | 246.9 ms | 287.2 ms | 0.8× | 0.86× | `NO-PATH` |
| `count` | 192.0 ms | 240.8 ms | 275.0 ms | 0.8× | 0.88× | `NO-PATH` |

`DISAGREE`, `SKEW`, `ND-SLOW` and `ND-UNSUPPORTED` are all zero on this set.

`MemberQ`/`Position`/`Count` stay `NO-PATH` and that is the honest reading: they
now give the *right* answer on every surface, and they still have no buffer fast
path — they materialise, which is register C.7, a performance gap rather than a
correctness one. Being slightly slower than the plain-List column is the cost of
the materialise; closing it means a buffer-aware scan for the literal-pattern
case, which is a separate piece of work.

- Full suite: 396 binaries, only the two documented pre-existing failures
  (`simplify_tests`' `Simplify[(Sqrt[x^2] - 1/Sqrt[x^2])/x^2]`, which fails
  identically on unmodified `HEAD`, and `crc_corpus_tests`' known hang).
- `make check-c99` and `make check-packed-aware` both pass; the latter with one
  fewer `EXEMPT` and one fewer `NOT_AWARE` entry than before.
- `Mod`/`Quotient` values checked against `wolframscript` element for element,
  including both sign conventions and the `Real`-in/`Integer`-out narrowing.
