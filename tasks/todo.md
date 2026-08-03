# Task: Implement `Ordering` (2026-08-03)

`Ordering[list]` → the 1-based permutation that sorts `list`
(`list[[Ordering[list]]] === Sort[list]`). Forms: `Ordering[list, n]`/`-n`
(n smallest/largest), `Ordering[list, seq]` == `Take[Ordering[list], seq]`,
`Ordering[list, seq, p]` (order by `p`). Any head; Association orders by values.
Required: highly efficient (packed fast path, consistent with `Sort`), and
supported in `Compile[]`/auto-compilation.

## Plan / progress

- [x] `builtin_ordering` in `src/sort.c` (+ `src/sort.h`): borrowed-subject
      stable index argsort (index tie-break over `qsort` → ties by original
      position), Take-spec via the shared `get_seq_spec_indices`, Association by
      value, any head → always a `List` of positions, `pack_offer`ed.
- [x] Factored the `p[a,b]` comparison into shared `cmp_with_p` (used by
      `custom_compare` and the ordering comparator).
- [x] Exported `get_seq_spec_indices` from `src/list/take_drop.{c,h}`.
- [x] Packed fast path: `ndstruct_ordering` (`src/ndstruct.{c,h}`) +
      `nd_argsort_i64`/`nd_argsort_real` stable merge argsort
      (`src/ndreduce.c`, `src/ndarray_internal.h`). int64-exact past 2^53.
- [x] `pack.c`: `Ordering` added to `AWARE` and `INT64_OK`.
- [x] Compile: `ND_FNS` row `{ "Ordering", ndstruct_ordering, 0, 0, true }`;
      new `NdFnSpec.int_result` → `nd_fn_result` forces `CT_INT` (a permutation
      is integer for any input dtype). Auto-compile free.
- [x] Registered `Protected` in `src/core.c`; `SYM_Ordering` in
      `src/sym_names.{c,h}`; docstring in `src/info.c`.
- [x] Docs: `## Ordering` in `docs/spec/builtins/structural-manipulation.md`;
      changelog entry in `docs/spec/changelog/2026-08-03.md`.
- [x] Tests: `test_ordering` (`tests/test_sort.c`),
      `test_ordering_builtin_matches_plain` (`tests/test_packed_list.c`),
      `test_cf_ordering` (`tests/test_compiledfunction.c`). All pass.

## Review

- Clean `-std=c99` build; `make check-c99` passes. No new POSIX/`long long`
  usage; all buffers `int64_t`.
- Every spec example reproduced verbatim in the REPL, including
  `Ordering[{2,6,1,9,1,2,3}, All, Greater]` → `{4,2,7,6,1,5,3}`.
- Audits: `check-c99`, `check-packed-aware` (Ordering in AWARE/INT64_OK),
  `check-array-exactness` (0 MIXED), `check-nd-surfaces`, and
  `check-compile-coverage` (Ordering compiles, no new gaps) all pass.
- `check-fastpath-sweep` exits 1 on a **pre-existing** stale `OFF_BUFFER`
  baseline (65 unrelated heads incl. `Order`, added earlier today, `D`,
  `Variables`, `TrueQ`, …). `Ordering` is NOT among them — `MATHILDA_PACK_DIAG=gate`
  shows it materialises **0** packed arguments (buffer path taken). My only
  `pack.c` change adds Ordering to AWARE, which can only reduce materialisation,
  so the red gate is not a regression from this task. (Recording the backlog in
  `OFF_BUFFER` / giving those heads buffer paths is a separate cleanup.)
- valgrind: leak numbers byte-identical to a no-Ordering baseline → zero leaks.
- Note: for a custom **strict** comparator (e.g. `Greater`) with tied elements,
  tie order is `qsort`-defined exactly as in `Sort[list, p]`; the guaranteed
  contract is `list[[ord]] === Sort[list, p]`. The default order is fully stable.

---

# Task: Implement `Order[e1, e2]` (2026-08-03)

`Order[e1, e2]` → `1` if e1 is before e2 in canonical order, `-1` if after, `0`
if identical. Structural (same order as `Sort`), not numerical:
`Order[6, Pi] == 1`, `Order[6, N[Pi]] == -1`. Also required: compilable inside
`Compile[]` and auto-compilable.

## Plan / progress

- [x] Interpreter builtin `builtin_order` in `src/sort.c` — one `expr_compare`
      call, sign inverted (`expr_compare < 0` ⇒ before ⇒ `+1`); shares the
      comparator every sorting routine uses. Arity ≠ 2 → `builtin_arg_error`.
- [x] Declared in `src/sort.h`; registered `Protected` in `src/core.c`.
- [x] `SYM_Order` in `src/sym_names.{c,h}`.
- [x] Docstring in `src/info.c` (terse, no examples).
- [x] Compile: `Order[a,b] = Sign[b - a]` via existing `SUB`+`SIGN` opcodes;
      branches added to both `infer_type` and `emit_node` in
      `src/compile/compile.c`, result type `CT_INT`. No new opcode. Complex/
      array/bool args decline → interpreter.
- [x] Auto-compilation: free (autocompile.c uses the same `compile_expr_ex`).
- [x] Tests: `test_order` in `tests/test_sort.c`; `test_cf_order` in
      `tests/test_compiledfunction.c` (incl. `CompileDiagnostics` Compiled→True
      and `Head[...]==Integer` parity).
- [x] Docs: `## Order` in `docs/spec/builtins/structural-manipulation.md`;
      changelog entry in `docs/spec/changelog/2026-08-03.md`.

## Review / results

- Build clean (`make`), `make check-c99` clean, `make check-compile-coverage`
  green ("no new gaps" — Order has no fast path so it's correctly outside the
  audit inventory).
- REPL: `{Order[a,a],Order[a,b],Order[b,a]} == {0,1,-1}`;
  `{Order[6,Pi],Order[6,N[Pi]]} == {1,-1}`;
  `Order @@@ Tuples[{0,1,2},2] == {0,1,1,-1,0,1,-1,-1,0}`; `Attributes[Order] ==
  {Protected}`; arity errors emit `Order::argrx` and stay unevaluated.
- Compiled: real/int/mixed all correct; `CompileDiagnostics` → `Compiled→True,
  ResultType→Integer` (3 instructions, not bailed); `Head` is `Integer`;
  `Table[Order[i,3],{i,1,5}] == {1,1,0,-1,-1}`.
- `sort_tests`, `compiledfunction_tests`, `list_tests` all pass.
- valgrind: no leak stack implicates `builtin_order`/`expr_compare`/`test_order`
  (only macOS objc/dyld framework baseline noise remains). NB: the literal
  `Order @@@ Tuples[...]` surfaced a **pre-existing** leak in `builtin_tuples`
  (unrelated to Order); the leak-sensitive test uses the explicit tuple list
  instead, and the literal Tuples form is covered by the REPL smoke test.

---

# Task: the sixth sweep — every builtin, packed arrays and `Compile[]` (2026-08-02)

`Commonest` was found costing 880 ms where `Tally` of the identical 10⁷ buffer
cost 21.5 ms, **after five sweeps and with four audits green**. That is the real
finding: the audits share a blind spot, and it is structural, not an oversight.
Two of them read the *source* (so a head with no dispatch is invisible), one
reads *element heads* (Commonest's were right), and the measuring one works from
a hand-written probe list with no `commonest` entry.

Second requirement, added mid-task: every builtin with a numeric fast path
should also be compatible with `Compile[]` and auto-compilation.

## Plan

- [x] **Build the audit that names nobody.** `tools/nd_fastpath_sweep.py`:
      enumerate every live symbol (`Names["*"]`, not the C scan — 313 heads are
      defined in `src/internal/*.m`), discover by TRIAL which of 12 call shapes
      each accepts, then either count materialisations under
      `MATHILDA_PACK_DIAG=gate` (the detector: a count, not a duration) or time
      them packed and under `MATHILDA_NO_PACK=1` (the severity).
- [x] Harden it against its own hazards: dimension-spec arguments
      (`PadLeft[vi, wi]` asks for `prod(wi)` elements and hung the first run),
      rendering heads by suffix rather than by name, tight timeouts.
- [x] **Build the Compile-coverage audit.** `tools/compile_coverage.py`: join
      the kernel registry and `pack.c`'s `AWARE` against the binary's own
      `CompileDiagnostics` over every typed argument shape.
- [x] Make its probe list honest — functional spellings (`H[f, v]`) and integer
      arrays, or `Map`/`Select`/`Fold`/`GCD` report the PROBE's gap as theirs.
- [x] Fix the Compile gaps that fall into coherent classes:
  - [x] array → scalar reductions (`OP_V_NDRED`): `Mean`, `Median`, `Variance`,
        `StandardDeviation`, `RootMeanSquare`, `Max`, `Min`
  - [x] array → array delegations (`ND_FNS`): `Differences`, `Ratios`, `Most`,
        `Rest`, `Clip`, `RotateLeft`, `RotateRight`, `MovingAverage`,
        `MovingMedian`, `TakeLargest`, `TakeSmallest`
  - [x] narrowing kernels over an array: `IntegerPart`, `UnitStep`
  - [x] exact-integer binary kernels over an array: `Mod`, `Quotient`, `GCD`,
        `LCM`, `ArcTan[v, y]`, `DivisorSigma`
  - [x] integer-only unary kernels over an array: `MoebiusMu`, `EulerPhi`,
        `IntegerLength`
- [x] Fix the two WRONG ANSWERS the audit turned up (below).
- [x] Ratchet `check-compile-coverage` against a checked-in `BASELINE`, so a new
      gap fails while the standing 52 are reported.
- [x] Tests, docs, changelog, `make` targets.
- [x] Fix `Subtract` / `Divide` — one rewrite away from the buffer — and the two
      wrong answers that surfaced while verifying them.
- [x] Fix `Range`'s silent truncation at 10^6 elements (raised mid-task by the
      user after the sweep tripped over it).
- [x] Run the gate pass over every builtin, verify each finding one-per-process,
      and record the 48 that survive in the tool's `OFF_BUFFER` ratchet.
- [ ] Work that queue. `Position` and `Count` first — the two a numeric workload
      reaches most, and both a single pass over the buffer.

## Review

### Two wrong answers, found by the Compile audit

**`Max[v]` lowered to the identity.** The scalar pairwise fold accepted an array
operand unchecked, and an empty fold returns its accumulator:

    Compile[{{v, _Real, 1}}, Max[v] + 1.][{3., 1., 7., 2.}]
      was  {4.0, 2.0, 8.0, 3.0}          interpreter: 8.

`Max[v]` alone was rejected downstream, so the wrong answer needed the head
inside a larger expression before it showed.

**The narrowing kernels declared the wrong element type.** `Floor` over a Real
array declared Real and wrote an `NDT_INT64` buffer. The declared type is what
downstream opcodes read the slot as, so:

    Compile[{{v, _Real, 1}}, Total[Floor[v]]][{1.5, 2.5, 3.5}]
      was  2.96439*10^-323               interpreter: 6   (that is the int64 6)
    Compile[{{v, _Real, 1}}, Floor[v] + 1][{1.5, 2.5, 3.5}]
      was  {2.0, 3.0, 4.0}               interpreter: {2, 3, 4}

Both are pinned in `tests/test_compiledfunction.c`.

### Two stale test expectations, fixed alongside

`PACK_MIN_ELEMENTS` moved 250 → 4 on 2026-08-02 and two assertions still named
the old boundary. The threshold test now *pins* it with
`pack_set_min_elements(250)` rather than inheriting it — a test about what
happens below a boundary must name the boundary.

### The packed surface

`Subtract` and `Divide` are now on `AWARE`: neither reads an element, both
rewrite to `Plus`/`Times`/`Power`, and the gate was firing at the head doing the
rewriting. Measured against a pristine build:

| | before | after |
|---|---:|---:|
| `Subtract[v, k]`, 10^6 float64 | 857.8 ms | 1.09 ms |
| `Divide[v, k]`, 10^6 float64 | 380.1 ms | 0.97 ms |
| `Subtract[iv, k]`, 10^6 int64 | 833.8 ms | 1.00 ms |

The infix `-` and `/` were already fast — the parser desugars them to `Plus` and
`Times`, so only the explicit heads ever reached the slow path. That is worth
recording because the gate diagnostic named `Subtract`, and the obvious reading
of that ("`v - 1.` is slow") was wrong.

Two wrong answers came out of verifying it. `Divide[packedList, 2.]` answered
`0.` — the Real branch fires when *either* operand is Real and then reads the
other through scalar type tests an array fails, defaulting to `0.0`. And
`packedList + x` was left **unevaluated** where the same short list threads;
that one was pre-existing and in `Plus`/`Times`/`Power`, but marking the two new
heads aware would have propagated it, so it had to be fixed at the source.

### The 48 heads still off the buffer

Verified one expression per head in a single process (the batched pass
over-attributes — see below) and recorded in `OFF_BUFFER`, grouped by what
closing them needs: element-wise and cheap (`Chop`, `Im`, `Rationalize`,
`Numerator`, `Denominator`, `MantissaExponent`, `RealDigits`); integer
predicates over an int64 buffer (`EvenQ`, `OddQ`, `CoprimeQ`, `PrimeNu`,
`PrimeOmega`, `LiouvilleLambda`, `JacobiSymbol`, `DigitSum`,
`IntegerExponent`); search and selection (`Position`, `Count`, `MemberQ`,
`Cases`, `FirstCase`, `Delete`, `ReplacePart`, `Level`, `Apply`, `Thread`);
grouping and ordering (`Split`, `Gather`, `GatherBy`, `PositionIndex`,
`Lookup`, `OrderedQ`, `MinimalBy`, `ReverseSortBy`); shape predicates
(`SquareMatrixQ`, `SymmetricMatrixQ`, `DiagonalMatrixQ`,
`UpperTriangularMatrixQ`, `ArrayFlatten`); comparisons, input side only
(`Less`, `LessEqual`, `Unequal`, `SameQ`, `UnsameQ`); and four that need an
array-against-array binary map before they can be done at all (`BesselJ`,
`BesselI`, `BesselK`, `QuotientRemainder`).

### Three ways the tool lied before it was believed

- `Head[result] =!= H` reads a `Listable` head's threaded list of *unevaluated*
  calls as a success. The test is now `FreeQ[result, H]`.
- `symtab_add_builtin` is not the symbol table — 313 live names have no C
  registration. The inventory is now `Names["*"]`.
- The gate counter is per PROCESS, and the pass runs 16 probes per process, so
  a head could be credited with a neighbour's materialisation (`Exp`, `Gamma`,
  `Variance`, `Clip` were all named wrongly). A verify phase re-asks every
  survivor one-per-process.

### Where the Compile surface stands

80 of the 234 heads inventoried had a numeric fast path and no compiled
lowering; **52 remain** (of 236 — `Subtract` and `Divide` joined the inventory
when they became packed-aware), grouped in `BASELINE` by what closing them
needs: an array×array
opcode (linear algebra, transforms), `ND_FNS` entries whose extra argument is
not an integer (`Join`, `Riffle`, `Partition`, `Pad*`), matrix producers, a
callback lowering (`Scan`, `Map*`, `NestWhile*`), and the SCALAR form of the
integer-only kernels. `Positive`/`Negative`/`NonNegative`/`NonPositive` are
exempt rather than open: they answer with a boolean List and there is no boolean
array dtype (§13 gap C.1).

Written up head by head in `COMPILE_MISSING.md`, which names the entry point
to delegate to where one exists. Two are not compiler gaps at all: `RowReduce`
and `NullSpace` are on `AWARE` and materialise the buffer on their own first
line, so they need a LAPACK path before a lowering means anything.

---

## Eigenvectors of an irrational algebraic eigenvalue (2026-08-03) — DONE

Closes the item left open by the 2026-08-02 SVD work: "a full-rank exact 3×3
returns *unevaluated* after 829 ms." That note guessed the failure was "in
`Eigenvectors` of the gram" — right location, wrong mechanism. `Eigenvectors`
did not fail to evaluate; it returned `{{0,0,0},{0,0,0},{0,0,0}}`, a **wrong
answer**, and the SVD then divided by a zero norm and declined.

- [x] Reproduce and localise: the bail is `sv_is_zero(nsq)` in step 5 of
      `svd_symbolic_core`, because the eigenvector really is the zero vector.
- [x] Root cause: `RowReduce[m - λI]` returns the **identity** for a matrix
      whose `Det` is exactly 0. `is_zero_poly` is a polynomial-identity test
      and cannot see the algebraic relations among nested cube roots.
- [x] Fix: `eigen_null_space_algebraic` — a column of `adj(m - xI)` reduced
      mod the minimal polynomial. No zero test needed during the computation.
- [x] Avoid the 15.9 s/eigenvalue `MinimalPolynomial` call: when the
      characteristic polynomial is irreducible it *is* the minimal polynomial
      (`IrreduciblePolynomialQ`, 0.16 ms). 47 s → 0.29 s.
- [x] Verified: SVD reconstructs to 1.7e-34, `UᵀU`/`VᵀV` orthonormal to 1e-33,
      singular values agree with LAPACK. Residuals ~1e-32 for the eigenvectors.
- [x] Edge cases held: defective zero-padding, identity/zero/rational/complex
      eigenvalues, symbolic matrices, repeated irrationals, generalised form.
- [x] 13 test suites green (rebuilt with the CORRECT target names — the first
      run used stale binaries); both new tests confirmed to fail without the fix.
- [x] Valgrind: leak summary byte-identical to the `Print[1]` baseline.
- [x] `make check-c99` clean. Docs + changelog updated.

**Still open, carried forward from 2026-08-02:**
- A packed SVD argument is 37% SLOWER than a plain one (56.7 vs 41.5 ms,
  300×300 float64) — it is on `AWARE` and then calls `linalg_delist_and_reeval`.
- A *reducible* characteristic polynomial with an irrational factor is now
  correct but pays `MinimalPolynomial` (≈33 s for a 4×4). Identifying which
  factor `λ` belongs to needs the zero test that is missing; solving it would
  mean threading the eigenvalue↔factor pairing out of `eigen_solve_poly`,
  which already factors before solving.
- A repeated irrational eigenvalue vanishes the whole adjugate (rank ≤ n-2);
  that eigenspace still gets zero-vector padding.
- `Diagonal` is still not implemented at all.
