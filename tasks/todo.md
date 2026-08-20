# Efficient Diophantine box search: Frye / Lander–Parkin / general fallback

Plan: `/Users/user/.claude/plans/we-should-have-a-noble-crescent.md`

## Phase 0 — Understand exact internals
- [ ] Read SICtx / SearchState structs, emit_full, si_verify, MitmEntry, build_result
- [ ] Read the full Booker template (si_solve_three_cubes_booker + helpers)
- [ ] Confirm ordering-constraint access (ord_a/ord_b, longest chain, effective_bounds)

## Phase 1 — Component 2: si_solve_separable_mitm (Lander–Parkin + separable)
- [ ] Shape-match separable additive single eq, fully bounded, engage-gate
- [ ] Ordering-aware chain split (hash prefix / iterate suffix), ordered generation
- [ ] __int128 sums; reuse MitmEntry sorted-array + binary search
- [ ] Register in si_try_special_forms; witness/decline contract
- [ ] Verify: Lander–Parkin returns tuple + k=2..6 multiples

## Phase 2 — Component 1: si_solve_biquadrate_frye (headline)
- [ ] Shape-match x^4+y^4+z^4==w^4 (unit coeffs, sign/perm normalise), all bounded
- [ ] Modular sieve tables: mod 625, mod 16, mod 9/13/29; 8k+1 factor constraint
- [ ] z-window (0.76w, w); residue-restricted decompose; int128 4th-root test
- [ ] find-first ascending; env-tunable node cap; witness/decline contract
- [ ] Verify: finds 95800/217519/414560/422481; measure & report timing

## Phase 3 — Component 3: si_solve_box_modsieve (non-separable)
- [ ] Modular-pruned nested enumeration; auto modulus M from exponents
- [ ] Register; witness/decline contract
- [ ] Verify vs Python brute force on boxed cross-term case

## Phase 4 — Tests, docs, gates
- [ ] tests/test_solve_integers.c new cases; regression suite green
- [ ] benchmarks/87 heldout.py cases (OK / proven {} / DECLINE, zero WRONG)
- [ ] docs/spec/builtins/solutions-of-equations.md + changelog 2026-08-17 + SOLVE_INTEGERS.md
- [ ] make check-c99; valgrind delta 0; make check-diophantine-heldout

## Review

All three components implemented in `src/solveint.c`, registered in
`si_try_special_forms` (Frye, separable MITM) and at the leaf-search decline
point (modsieve).

**Results (measured on this machine):**
- **Lander–Parkin** `x^5+y^5+z^5+w^5==r^5 && 0<x<y<z<w<r<1000` → complete set
  (144-tuple + its 2×–6× multiples), **~6.5 s**. Ordering-aware, `__int128`,
  modular residue sieve. (Was: unevaluated.)
- **Frye** `x^4+y^4+z^4==w^4 && 0<x<y<z<w<10^6` → `{95800,217519,414560,422481}`,
  **~11 min single-threaded** (654 s, 62.6e9 decompose iters). Full modular +
  factor sieve. (Was: unevaluated. Frye 1988 needed a supercomputer.)
- **General non-separable** big boxes (e.g. `x^2+xy+y^2==z^2 && …<15000`) →
  complete set via modular-sieved exhaustive leaf search, **~6 s**. (Was:
  unevaluated.) Separable general boxes (`x^2+3y^4==6z^4`) handled by the
  existing/MITM path.

**Verification:** 3 new unit tests pass; held-out gate `three-cubes-eq-cube-30`
OK, 0 wrong answers; modsieve proven identical to the ordinary engine on the
overlap sub-box; macOS `leaks` 0 bytes; `make check-c99` clean.

**Key design decisions (with the user):** witness semantics (find-and-return,
never an unproven `{}`); single-thread + factor sieve for Frye; general fallback
covers non-separable. All engage-gated so existing small-box behaviour is
unchanged (regression-safe).
