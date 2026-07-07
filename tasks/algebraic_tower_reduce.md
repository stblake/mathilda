# Task: rigorous FLINT reduction over genuine-algebraic parametric towers

## Problem
`Simplify[D[F,x] - integrand]` (F = Goursat cube-root integral) does not reduce to
0. The expression is a rational function over the tower
`K = Q(x,k)(α, β, ω)` with
- α = (x(1-x)(1-kx))^(1/3),  α³ = x(1-x)(1-kx)   (genuine algebraic, poly radicand)
- β = k^(1/3),               β³ = k              (genuine algebraic, symbol radicand)
- ω = (-1)^(1/3),            ω²  = ω - 1         (root of unity, Φ_6)

Current FLINT bridge handles only `Sqrt[symbol]` (q=2, free-symbol radicand) via a
transcendental collapse. Cube-roots and polynomial radicands are rejected by
`ps_detect` and `extension_autodetect`, so this class falls to the non-FLINT
Phase-E / radrat paths, which under-reduce or blow up (OOM).

## Constraints (from user)
- Rigorous algebraic reduction only. **Simplify must never call PossibleZeroQ.**
- Leverage FLINT. Route Cancel/Together/Simplify through it.

## Design — Gröbner-ideal reduction over Q(params)[gens]
Represent K as `Q[params, gens] / I`, `I = (g_l^{q_l} - radicand_l, ...)`, monic in
each generator. With generators as the *leading* (first) LEX variables, the minpolys
are a Gröbner basis (distinct monic leading monomials → Buchberger 1st criterion), so
`fmpz_mpoly_divrem_ideal` gives the canonical normal form. VALIDATED via probe.

Zero test / Together: combine E → N/D over Q[params,gens] (via existing
`expr_to_mpolyq` after substituting each radical/root-of-unity with a fresh gen
symbol), reduce N mod I. `N ≡ 0` ⟺ `E ≡ 0` in K (sound always; complete when each
minpoly irreducible — our case). No inversion, no primes/CRT for the zero reduction.

## Milestone A — zero-reduction + Together numerator reduction (this session)
- [ ] A1. `genalg_detect`: collect genuine-algebraic gens (Power[base,p/q], q≥2,
      base poly/symbol — canonicalise base by Expand; (-1)^(p/q) → root-of-unity gen
      via Φ_{2q}); collect params. NULL if none / out of scope.
- [ ] A2. Generator substitution E → Ē (radical^p → gen^p, (-1)^(p/q) → ω^p).
- [ ] A3. `flint_algebraic_field_normalize(e)`: mpoly ctx (gens first),
      `expr_to_mpolyq(Ē)` → N/D, reduce N (and D) mod I via divrem_ideal. N→0 ⟹
      return 0; else render reduced N/D back and eval.
- [ ] A4. Wire into `flint_cancel_fraction` (rat.c) first-try on whole expression.
- [ ] A5. Tests + Goursat `D[F,x]-integrand // Simplify == 0`; regressions; valgrind.
- [ ] A6. Docs + changelog + plan-status update.

## Milestone B — full field arithmetic [DONE 2026-07-04]
User chose (AskUserQuestion) "GCD-cancel in Cancel + RootReduce": deliver BOTH the
WL-faithful GCD-cancel/Together over cube+ roots AND a rationalising RootReduce.

- [x] B1. `flint_algebraic_field_canonical(e)` — rationalised canonical form via
      field inversion. K = Q(params)-vector space (monomial basis over the gens);
      multiply-by-denominator is a Q(params)-linear map (products reduced mod I via
      `fmpz_mpoly_divrem_ideal`); solve `M_denR x = coords(numR)` over Q(params)
      (`gr_mat_nonsingular_solve` on `gr_ctx_init_fmpz_mpoly_q`). Denominator
      rationalised to Norm = det(M). Coordinates recombined over lcm of their
      denominators → single fraction. Validated first with a standalone FLINT probe.
- [x] B2. `RootReduce` builtin in its OWN file `src/rootreduce.c` (+ `.h`),
      registered in `core_init` after `rat_init`; Listable|Protected; docstring.
      Identity when no algebraic generator.
- [x] B3. `flint_algebraic_field_together(e)` — WL-faithful combine of a Plus of
      cube+ root fractions treating radicals as FREE kernels (no rationalisation),
      via `fmpz_mpoly_q`. Gated to genuine (symbol-radicand) gens, rou allowed,
      constant radicals deferred. Wired into `flint_cancel_fraction` gated to Plus
      (so it never pre-empts the relation-aware single-fraction GCD path).
- [x] B4. Tests `tests/test_rootreduce.c` (`rootreduce_tests`, added to CMake
      COMMON_SRC + as an executable): exact forms, value preservation via the
      independent ideal-reduction zero test agreeing with the matrix inversion,
      rationalised-denominator checks, idempotence, Together/Cancel, stress batch.
- [x] B5. Docs (algebra.md RootReduce + cube-root Together note), changelog. valgrind
      clean (0 stacks in new code; 13,440B/420 blocks are macOS dyld/libobjc noise).

### Review — DONE (Milestone B), 2026-07-04
`RootReduce[1/(1+k^(1/3))]` → `(1-k^(1/3)+k^(2/3))/(1+k)`;
`RootReduce[1/(a+b Sqrt[2]+c Sqrt[3])]` → norm-denominator canonical form;
`Together[1/(x-k^(1/3))+1/(x+k^(1/3))]` → `2x/(x^2-k^(2/3))` (was left uncombined).
Key design: rationalisation (RootReduce) uses the ideal reduction; Together/Cancel
of sums uses FREE-kernel combination (WL-faithful, keeps radicals in denominator).
The Plus-gate keeps single-fraction relation-aware Cancel (`(x^3-k)/(x-k^(1/3))`)
intact. RootReduce fast; verification via full Simplify is slow on multi-radical
number fields, so tests use the cheap `RootReduce[RootReduce[e]-e]==0` self-check
plus the independent `flint_algebraic_field_normalize` cross-check on genuine gens.
No regressions: flint_bridge/rat/simp/Goursat/extension/radical/trigrat/intrat green;
simplify_tests 4 pre-existing FAILs unchanged; intrat_corpus 5 DIFF-NONZERO == baseline.

Milestone C (not scoped): rationalise inside Cancel/Together by default was
explicitly declined by the user (diverges from WL); RootReduce is the rationalising
entry point instead.

## Review — DONE (Milestone A), 2026-07-04
Implemented `flint_algebraic_field_normalize` (src/poly/flint_bridge.c, ~300 lines)
+ header decl + wired into `flint_cancel_fraction` (src/rat.c). All A1–A6 done.

Result: `Simplify[D[F,x]-integrand]`, `Together[…]`, `Cancel[…]` → 0 rigorously
(FLINT Gröbner-ideal reduction, no PossibleZeroQ). Soundness verified: nonzero /
perturbed / `Sqrt[x^2]-x` cases stay unreduced (never wrongly 0); the `g²=x²`
relation cannot fake `Sqrt[x^2]=x`.

Two bugs caught + fixed during dev:
1. Params that appear ONLY inside a radicand (e.g. k in k^(1/3)) vanish from the
   substituted `ebar` but the relation `gen^Q - radicand` still needs them → must
   collect symbols from each generator's `base_exp` into `fvars`, not just `ebar`.
2. Test global-symbol pollution: `d = D[F,x]-f` clobbered the bare `d` a later
   `FLINT`Det[{{a,b},{c,d}}]` test read → use unique names (gaF/gaf/gad).
Plus a 0/0 soundness guard (return 0 only if numerator≡0 AND denominator≢0 mod I).

Regressions (all vs baseline): rat/intrat/radical_*/simp/simp_algebraic_cuberoot/
integrate_goursat/extension_*/fullsimplify/trigrat green; flint_bridge +2 tests
green; simplify_tests 4 pre-existing FAILs (nested-radical denesting, unchanged,
return nonzero not 0); intrat_corpus 5 DIFF-NONZERO == HEAD baseline (proved by
building HEAD and getting the identical 5). valgrind: 0 leak stacks in my code
(the 64KB/626 blocks are pre-existing classical together_recursive/poly_gcd leaks).
