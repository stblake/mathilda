---
title: qa substrate — accept Power[c, p/q] for general integer p
date_started: 2026-05-13
status: shipped 2026-05-13 (Phases A–E); deeper Together[D[Integrate[...]]] case deferred
---

## Status (2026-05-13)

Phases A (recogniser broadening), B (lift inputs), C (complexity gate),
D (Q(γ)-arithmetic in the no-variable path), and E (input-side radicand
canonicalisation) all shipped with zero regressions
(96 pass / 14 fail — identical pre-existing failure set).

Key user-visible wins:
- `Cancel[Power[2,-2/3] x^2 - x/Power[2,1/3], Extension -> Power[2,1/3]]`
  → `1/2 (-2^(2/3) x + 2^(1/3) x^2)` (Phase B lift)
- `Cancel[(Sqrt[2]+Sqrt[3])(Sqrt[2]-Sqrt[3]), Extension -> Automatic]`
  → `-1` (Phase D Q(γ)-arithmetic, replaces the degree-16 γ-polynomial
  expansion that the previous code produced)
- `Together[Sqrt[-1/9/2^(2/3) + 2/9·2^(1/3)] - Sqrt[1/3/2^(2/3)],
  Extension -> Automatic]` → `0` (Phase E input-side canonicaliser
  collapses mathematically-equal radicands before tower substitution).
  Test: `tests/test_simp_algebraic_cuberoot.c::test_together_equal_nested_radicals`.

Remaining work:
- The deeper motivating example
  `Together[D[Integrate[a x/(x^3+2), x], x], Extension -> Automatic]`
  still produces a non-collapsed form containing `Sqrt[1/3/2^(2/3)]`
  embedded in arithmetic.  Closing this requires recognising
  `Sqrt[1/3/2^(2/3)] = Sqrt[3]/(3·2^(1/3))` (factoring the integer
  power out of the radicand) so the autodetected generator set
  contains both `Sqrt[3]` and `2^(1/3)` rather than the opaque
  `Sqrt[1/3/2^(2/3)]`.



# Plan

Extends picocas's algebraic-number substrate to accept `Power[c, p/q]` with
general integer `p` (currently only `p == 1` is accepted by G8's
`expr_is_atomic_algebraic`).  Closes the motivating example
`Together[D[Integrate[a x/(x^3+2), x], x], Extension -> Automatic]` →
`a x^2/(x^3+2)`.

## Design decisions (locked)

- **Q1**: A — `qa_resolve_extension` returns `*render_out` in the
  natural form `Power[c, 1/q_red]` for the new general-p case.  When
  `p == 1` (existing behaviour), `*render_out` is bitwise-identical
  to what is returned today (preserves `Sqrt[c]` and `Power[c, 1/n]`
  shapes byte-for-byte).
- **Q2**: yes — accept both positive and negative `p` uniformly via
  gcd reduction.
- **Q3**: falls out of Q1.A automatically (auto-detect merges via
  expr_eq on the natural form).
- **Q4**: reject `c ∈ {-1, 0, 1}` in the new branch; non-integer base
  + `p != 1` falls through to existing behaviour (no G8 rewrite).

## Phases

### A. Recogniser broadening
- `src/qafactor.c::expr_is_atomic_algebraic`: extend Power branch to
  accept any `Power[c, p/q]` with `q_red >= 2` after gcd reduction.
- `src/qafactor.c::qa_resolve_extension`: in Power branch, accept
  general p; return QAExt for `y^q_red - c` with `*render_out` set to
  the natural form per Q1.A.  Keep `p == 1` shortcut for byte-for-byte
  identical existing-test output.

### B. Lift inputs containing `Power[c, p/q]`
- **B1** `qa_alpha_power_signed(ext, p)` (qa.{c,h}): α^p for any
  integer p, using `qa_inverse` for p < 0; repeated-squaring for
  efficiency.
- **B2** `expand_radicals_to_atomic_poly(poly, c_base, q_natural,
  alpha_sym, ext)`: walk poly; rewrite every `Power[c_base, p/q]`
  whose reduced denominator divides q_natural into the polynomial-
  in-alpha_sym form `qa_alpha_power_signed(ext, p_red · q_natural/q_red)`.
  Raw expr substitution — must NOT call `evaluate` (would re-canonicalise
  via the Times absorber).
- **B3** Wire into `qa_expr_to_qaupoly_with_alpha`: add `c_base` and
  `q_natural` parameters (option B3.a from the plan).  Lift callers
  pass `c_base = 0` (sentinel) for non-radical-extension cases.
- **B4** G6 tower path: stash `(c, q)` per generator in `QATower`;
  apply preprocessing per-generator inside the substitute-α_i loop
  in `qa_factor_with_extension_tower` and `qa_cancel_with_tower`.

### C. Autodetect verification (no code change anticipated)
- Verify autodetect_walk's existing `q >= 2` accept still works.
- Verify autodetect_build_tower → qa_resolve_extension → Phase A path.
- End-to-end smoke test on the motivating example.

### D. Polish
- D1: keep `α^k with k >= ext->deg` safety belt in lift.
- D2: refuse non-integer base + `p != 1` in qa_resolve_extension's G8
  fallthrough (documented limitation).
- D3: `qa_factor_inner` lift in the G5 path needs B preprocessing too.

## Safety-check corpus (must keep working unchanged)

```c
/* 1. Pure G5, p == 1, Sqrt: most-trodden path. */
Factor[x^2 - 2, Extension -> Sqrt[2]]
    →  (x - Sqrt[2]) (x + Sqrt[2])

/* 2. G5, p == 1, cube root: user-form Power[2, 1/3]. */
Factor[x^3 - 2, Extension -> Power[2, 1/3]]
    →  (x - 2^(1/3)) (x^2 + 2^(1/3) x + 2^(2/3))

/* 3. Together with autodetect — the headline. */
Together[D[Integrate[a x/(x^3+2), x], x], Extension -> Automatic]
    →  a x^2 / (x^3 + 2)

/* 4. G6 tower with Sqrt[1 + Sqrt[2]] autodetect → degree 4. */
extension_autodetect("Sqrt[1 + Sqrt[2]]") →  n=1, deg=4

/* 5. Cancel with explicit Power[2, p/q] in input — Phase B path. */
Cancel[(x^3 - 2) / (x - Power[2, 1/3]), Extension -> Power[2, 1/3]]
    →  x^2 + 2^(1/3) x + 2^(2/3)
```

Plus full existing-test invariance: 94 PASS / 16 FAIL must remain.
