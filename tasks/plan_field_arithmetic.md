# Plan: native number-field arithmetic for ParallelMixedTower — solve the missing Charlwood integrals at CAS speed

## Goal (from the user, verbatim intent)

Solve the ParallelMixedTower "missing" integrals **at least as fast as other CAS**,
with **robust algorithmic fixes/extensions — no hacks** (no invocation gates, no
"explicit-vs-cascade" tricks). A fast algebraic-field path is correct on *every*
call site, so it *also* dissolves the DSolve Q(i) regression automatically (a fast
Q(i) assembly cannot blow up).

## Root cause (profiled, not guessed)

The ansatz linear system is assembled over one number field K = Q(θ) with
`AlgebraicNumber[θ, coords]` coefficients, but **all the field arithmetic runs
through the WL evaluator on the expression representation**, with an expr↔qqbar
round-trip per scalar operation. `sample` of DSolve[1219] (Q(i), the worst case)
over an 8 s window, max self-samples per function:

| function | samples | what it is |
|---|---|---|
| `expr_expand_impl` | 592 | `Expand` of field-coefficient polynomials |
| `multiply_two` | 577 | `Times` inner loop over field terms |
| `flint_qqbar_algebraic_number` | 354 | **re-canonicalises** every `AlgebraicNumber[…]` the evaluator re-evaluates (`to_qqbar` + `algint_generator` + `poly_to_algnum`) |
| `poly_to_algnum` | 330 | builds `AlgebraicNumber` from an `fmpq_poly` result |
| `qqbar_to_expr` | 282 | rebuilds θ's `Root`/`Sqrt` expr per result |
| `fmpq_mpoly_pow_ui` / `to_mpoly` / `mpoly_ctx_init` | ~200 | the FLINT rational-mpoly bridge (Together/Cancel/Expand) |

The dominant cost (≈ 40 % of samples) is `flint_qqbar_algebraic_number` +
`poly_to_algnum` + `qqbar_to_expr`: **per-operation canonicalisation of field
elements through expr and qqbar.** A16 (a genuine Charlwood radical) is *small* and
so solves over Q(i) in 0.84 s; DSolve[1219] is *large* (a complex-exponential
intermediate) and pays this per-op cost thousands of times → 40 s. **Q(i) is not
inherently slow — the per-op overhead is.** This is exactly what SymPy avoids with
its `QQ<alpha>` domain and Mathematica with compiled AlgebraicNumber arithmetic.

### Why point-fixes are not enough (measured)
- **mpoly context cache** (landed, `flint_bridge.c`): removes `mpoly_ctx_init`
  (#1 leaf); DSolve[1219] 40.9→34.7 s (~15 %). Real but partial.
- **qqbar_to_expr Root memo** (landed, degree ≥ 3): helps the *compositum* cases,
  but Q(i)'s θ = I is **degree 2**, so it does nothing for DSolve[1219] (0 % there).
- No single memo on the existing pipeline reaches parity: the arithmetic itself runs
  through `Expand`/`Times`/`Plus` on exprs, which re-evaluate (and re-canonicalise)
  every field subexpression.

## The robust fix — native K = Q(θ) coefficient arithmetic

Represent the assembly's field-coefficient polynomials **natively** and do the hot
arithmetic in C over the number field, converting expr↔native **once per rung**
(not per operation). FLINT antic `nf_elem` (already linked, `nf_init`/`nf_elem_*`)
is the field element; a polynomial in the generators over `nf_elem` is the object.

Phased so each phase is measurable and independently correct:

### Phase 0 — foundation (LANDED this session, keep)
mpoly context cache; `wl_root_index` per-minpoly memo; `qqbar_to_expr` degree≥3
memo; guarded FLINT allocator (honest budget); branch-resolution gate (G3).
All verified correct (poly/algebra suite green; RootReduce byte-identical).

### Phase 1 — a native "field element" fast path for the scalar arithmetic
Cut the ≈ 965-sample canonicalisation cost without changing the .m's structure:
1. **`nf` context cache** keyed by the minimal polynomial of θ (like the mpoly
   cache): `nf_init` once per field, not per op.
2. **`AlgebraicNumber` +/*/^ via `nf_elem`** (`plus.c`/`times.c`/`power.c` →
   `flint_qqbar_algnum_*`): operate on the coords as `nf_elem` mod the cached
   minpoly — no `to_qqbar` of θ, no `qqbar_to_expr` rebuild, no re-canonicalisation.
   The generator expr is threaded through unchanged (θ is fixed and canonical), so
   the result is `AlgebraicNumber[θ, coords']` with the *same* θ expr — identical
   output, no round-trip.
3. **Canonical fast path in `flint_qqbar_algebraic_number`**: if `coords` are
   already reduced (length ≤ deg minpoly, read from the cached `nf`), return the
   input unchanged — skip `to_qqbar`/`poly_to_algnum` on re-evaluation.
- **Gate:** dedicated canonical-arithmetic test (add/mul/pow/zero over Q(i),
  Q(√2), a Root field) byte-identical to today; then measure DSolve[1219] and A16.
  Expected: the 965-sample block collapses; DSolve[1219] toward its Q time (~5 s).

### Phase 2 — native field-coefficient *polynomial* arithmetic (if Phase 1 short)
The remaining cost is `Expand`/`Times` over polynomials in the generators with
`nf_elem` coefficients (592 + 577). Provide a C builtin the .m assembly calls for
its hot step — `ParallelMixed`Private`FieldExpand[poly, gens, θ]` (and a matching
`FieldCoefficientRules`) — that expands/reads a generator-polynomial with
`AlgebraicNumber[θ,…]` coefficients using a native `nf_elem`-coefficient multivariate
polynomial (FLINT `gr_poly`/`gr_mpoly` over the `nf` ring, or a compact custom
representation). Convert expr↔native at the call boundary only.
- **Gate:** `FieldExpand`/`FieldCoefficientRules` byte-identical to `Expand`/
  `CoefficientRules` on field-coefficient inputs (fuzz over random Q(θ) polys);
  then re-measure.

### Phase 3 — extensions for the still-unsolved cases
1. **P8** — `ToNumberField` over a **conjugate-root compositum** (the splitting
   field of the quartic special). Re-land the precision-escalating `in_field`
   (this session's G2) but **degree-gated** (escalate only for deg > 6) so it never
   slows the low-degree common path (the earlier ungated version regressed dsolve).
   With Phase 1/2, P8's field arithmetic is fast, so the compositum solve fits the
   budget.
2. **Degree-16 cases (P4, A2, A3, A27, A35, A40)** — genuine norm/S'-unit search
   over a degree-16 field. Re-profile after Phase 1/2; if still over budget,
   the remaining cost is the search itself — treat each as its own algorithmic item
   (e.g. reuse the field factorisation, prune the S'-unit ansatz). These are the
   conic ladders whose specials split over Q(√5, √(1±√5)) and Q(√3, i).

### Phase 4 — verification & release
- Full unit suite green **including the dsolve corpus** (the regression is gone
  because Q(i) is now fast — this is the acceptance test for "no blowups").
- Charlwood 50 at the paper's protocol: every solved case ≤ its Mathematica time
  where the field work dominates; no case slower than before.
- `make check-c99`, valgrind on the new C paths, changelog, `$VersionNumber` bump,
  tags.

## Design notes / correctness invariants
- **θ is fixed and canonical per field** (minted once by `ToNumberField`), so
  threading the generator expr through arithmetic is sound — no need to re-derive
  it. The `nf` context is a pure function of θ's minimal polynomial.
- **Exactness is preserved**: `nf_elem` arithmetic is exact field arithmetic mod
  the minimal polynomial — identical to today's `fmpq_poly` + `rem`, just without
  the expr/qqbar wrapping. The "a zero is exactly 0" invariant the assembly relies
  on (`CoefficientRules` dropping zero monomials) holds because `nf_elem` reduces
  to the canonical zero.
- **No behaviour change off the field path**: pure-rational and single-radical
  towers never build an `nf` context; the fast path is entered only when an
  `AlgebraicNumber` is present.

## Risk & sequencing
- Phase 1 is the highest-leverage, lowest-risk step (localised to the AlgebraicNumber
  arithmetic in `flint_qqbar.c` + `plus/times/power.c`, guarded by a byte-identical
  canonical-arithmetic test). Do it first, measure, and reassess whether Phase 2 is
  needed for parity.
- Phase 2 is a real C addition (native nf-coefficient polynomials); gate it behind a
  differential fuzz test vs `Expand`/`CoefficientRules`.
- Phase 3.1 (P8) is independent; Phase 3.2 (degree-16) may reveal that a couple of
  cases need algorithmic work beyond arithmetic speed — surface that honestly rather
  than force them under a budget.

## What is already landed (safe, regression-free) vs. gated on this plan
- **Landed & safe:** Phase 0 (mpoly cache, memos, honest budget, G3 branch fix).
  These only ever make cases faster/more-correct; never slower.
- **Gated on Phase 1 (the fast field arithmetic):** the `AlgAtoms` Complex detection
  that unlocks A1/A16/A19/A20/A37. It is correct but, until Q(i) is fast, it couples
  those wins to the DSolve regression. **Do not ship it before Phase 1 lands.**
