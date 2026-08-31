# DSolve — Hypergeometric/Kummer recognizer (§1c SpecialFunctionForm)

Add Kummer (₁F₁) and Gauss (₂F₁) recognizers to `DSolve`SpecialFunctionForm`.
Verifies for free via the existing HypergeometricPFQ derivative rule
(numerically confirmed: all three back-substitutions → 0.0).

## Implementation
- [ ] src/calculus/dsolve_specialform.c: Kummer block (P=b/x−1, Q=−a/x)
      → C[1] ₁F₁[a,b,x] + C[2] x^(1−b) ₁F₁[a−b+1,2−b,x]; decline if IntegerQ[b]
- [ ] src/calculus/dsolve_specialform.c: Gauss block (W=x(1−x), Q=−ab/W, P=(c−(a+b+1)x)/W)
      → C[1] ₂F₁[a,b,c,x] + C[2] x^(1−c) ₂F₁[a−c+1,b−c+1,2−c,x]; solve t²−St+Pr for a,b;
      decline if IntegerQ[c]
- [ ] Expand header comment + docstring (add verification "keep-on-undecidable" note)

## Tests
- [ ] tests/test_dsolve.c: t_method_hypergeometric_kummer, t_method_hypergeometric_gauss
      (pin method, Head===List, PossibleZeroQ residual); + integer-b decline test
- [ ] tests/test_dsolve_m5_stress.c: forward-generator family over non-integer a,b (,c)

## Docs
- [ ] docs/spec/builtins/calculus.md: extend SpecialFunctionForm row (~L676)
- [ ] docs/spec/changelog/2026-08-31.md: DSolve entry
- [ ] DSOLVE_PLAN.md: §1c SpecialFunctionForm note → Kummer/Gauss done

## Gates
- [ ] make -j build clean
- [ ] REPL spot-checks (Kummer, Gauss, pinned, ?docstring, integer-b decline)
- [ ] ctest -R dsolve
- [ ] valgrind leak-check (engine allocations balanced)
- [ ] make check-c99
- [ ] commit + push

## Review — DONE

All items complete. `DSolve\`SpecialFunctionForm` now recognises Kummer (₁F₁) and
Gauss (₂F₁) second-order linear equations (`src/calculus/dsolve_specialform.c`),
appended as two `if (!general)` blocks after Bessel; shared quadratic-root helper
`specialform_quad_roots` (FactorList linear factors → radical-free `a,b`, radical
`Solve` fallback) and `specialform_is_number` (the numeric-exponent gate).

**Verified:**
- Numeric Kummer/Gauss → closed form; symbolic-`a` Kummer & symbolic-`a,b` Gauss
  (numeric exponent) → clean closed form; Airy/Bessel regressions intact.
- `ctest -R dsolve` → 3/3 pass (unit + both stress suites); 6 new unit tests +
  2 forward-generator stress families (24 solves, all back-substitute to 0).
- `make check-c99` PASS.
- valgrind: recognizer + `quad_roots` add **zero** leaks (0 direct-alloc blocks,
  0 `quad_roots` frames); only pre-existing Simplify/engine baseline remains.

**Key finding — verification hang (worked around, not a regression):**
The second solution carries `x^(1−b)`/`x^(1−c)`. With a **symbolic** exponent the
verify residual is a symbolic-power + pFq sum on which `zero_test`/PossibleZeroQ
**hangs** (pre-existing). Gate: emit only when the exponent parameter (`b`/`c`) is
a `NumberQ` non-integer; symbolic/integer exponents decline to Frobenius. Confirmed
via `git stash` that fully-symbolic Gauss also hangs on the **pre-change** binary
(the hang is in Kovacic on symbolic coefficients — unrelated to this change).

**Follow-ups (out of scope):** affine/Möbius normalisation for non-canonical
singular points; a `zero_test` improvement (or verify time-box) to lift the
symbolic-exponent gate; LegendreP recognizer (needs LegendreP/Q derivative rules).
