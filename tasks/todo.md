# Fix "Binomial-class" numeric failures across Mathilda

## Group 1 — silently-wrong coercion (correctness bugs)
- [x] `Mod` — fixed (mod_operand_to_double + mod_set_mpfr); 2-arg/3-arg/MPFR all correct
- [x] `Divide` — fixed (guard machine path to machine-real operands; else exact rewrite)
- [x] `QuotientRemainder` — auto-fixed via `builtin_mod`; verified

## Group 2 — missing complex path
- [x] `Factorial` — complex arg → `Gamma[z+1]`

## Group 3 — no numeric path for non-integer args (existing functions)
- [x] `Factorial2`  = 2^(z/2) (Pi/2)^((Cos[z Pi]-1)/4) Gamma[z/2+1]
- [x] `FactorialPower[x,n]` = Gamma[x+1]/Gamma[x-n+1] (non-integer n)
- [x] `BarnesG` — asymptotic + recurrence as expression (real/complex/MPFR); coeff = B_{2k+2}/((2k)(2k+2))
- [x] `Hyperfactorial` = Gamma[z+1]^z / BarnesG[z+1]
- [n/a] `CatalanNumber`/`Multinomial`/`Subfactorial` — NOT implemented at all (missing features, out of scope)

## Cross-cutting
- [x] Thorough unit tests: tests/test_numeric_domain.c (12 groups, CMake-registered)
- [x] valgrind clean (identical to baseline); C99 gate passes; docstrings + spec + changelog

## Review
All Group 1 correctness bugs (Mod/Divide/QuotientRemainder) and Group 2/3 numeric
continuations (Factorial complex, Factorial2, FactorialPower, BarnesG, Hyperfactorial)
are implemented, each built as an expression reusing Gamma's complex+MPFR kernels.
BarnesG required a numerically-verified coefficient correction
(B_{2k+2}/((2k)(2k+2)), not the textbook (2k)(2k+1)(2k+2)). Exact/integer inputs
unchanged. CatalanNumber/Multinomial/Subfactorial excluded (not implemented at all).
Tests: numeric_domain_tests + core/bigint/modular all pass; valgrind clean; c99 gate green.
