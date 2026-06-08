# Task: Implement `Gamma` builtin in `src/gamma.c`

Faithful recreation of Mathematica's `Gamma`:
- `Gamma[z]` — Euler gamma function Γ(z)
- `Gamma[a, z]` — upper incomplete gamma Γ(a,z) = ∫_z^∞ t^{a-1} e^{-t} dt
- `Gamma[a, z0, z1]` — generalized = Γ(a,z0) − Γ(a,z1)

Attributes: Listable, NumericFunction, Protected.

## Design (dispatch in `builtin_gamma`)

### Gamma[z] (1 arg)
1. Exact integer / half-integer → reuse `Factorial[z-1]` machinery
   (Gamma[z] = (z-1)!): integers → (n-1)! exact/BigInt or ComplexInfinity
   for n≤0; half-integers → rational·Sqrt[Pi]. (Build `Factorial[z-1]`,
   eval; if it stays a `Factorial[...]` head, free + return NULL.)
2. Symbolic infinities: Infinity→Infinity, -Infinity→Indeterminate,
   ComplexInfinity→ComplexInfinity, Indeterminate→Indeterminate.
3. Machine real (EXPR_REAL) → `tgamma`. (poles → ComplexInfinity)
4. Machine complex (Complex[real,real], inexact) → Lanczos (double complex).
5. MPFR real → `mpfr_gamma` at input precision.
6. MPFR complex (Complex[mpfr,mpfr]) → Lanczos via exported `mpfr_complex_*`
   primitives (exp/log/sin + manual complex mul) with reflection for Re<1/2.
7. Otherwise NULL (stay symbolic).

### Gamma[a, z] (2 args)
1. Exact rewrites: Gamma[a,0]→Gamma[a]; Gamma[1,z]→E^-z; Gamma[a,Infinity]→0.
2. Numeric real a,z (machine or MPFR) → `mpfr_gamma_inc`
   (returns EXPR_REAL for machine inputs @53-bit, EXPR_MPFR for mpfr inputs).
3. Otherwise NULL.

### Gamma[a, z0, z1] (3 args)
Rewrite to `Gamma[a,z0] - Gamma[a,z1]` and evaluate.

## Files
- [ ] NEW `src/gamma.c` — builtin + helpers + `gamma_init()`
- [ ] NEW `src/gamma.h` — `void gamma_init(void);` + `Expr* builtin_gamma(Expr*);`
- [ ] `src/core.c` — `#include "gamma.h"` + call `gamma_init();`
- [ ] `src/info.c` — docstring for Gamma (terse, no examples)
- [ ] `src/sym_names.c` — already has `SYM_Gamma` (verify, no change)
- [ ] NEW `tests/test_gamma.c` — extensive coverage
- [ ] `tests/CMakeLists.txt` — add `gamma_tests` target
- [ ] `docs/spec/builtins/` — add Gamma entry (special-functions page)
- [ ] `docs/spec/changelog/2026-06-08.md` — changelog note (Mon of ISO week)

## Out of scope (future work)
- D[Gamma...] / Integrate[Gamma...] / Series[Gamma...].

## Verification
- [ ] `make -j` clean (`-std=c99 -Wall -Wextra`)
- [ ] Build + run only `gamma_tests` (scoped)
- [ ] valgrind diff vs Sin[1.0] baseline — no Mathilda-src leaks
- [ ] Spot-check REPL against spec example values

## Review

**Status: complete.** `Gamma[z]`, `Gamma[a,z]`, `Gamma[a,z0,z1]` implemented
in `src/gamma.c` (+ `gamma.h`), wired into `core.c`, docstring in `info.c`.
`SYM_Gamma` was already interned in `sym_names.c` (no change needed).

Files touched:
- NEW `src/gamma.c`, `src/gamma.h`
- `src/core.c` (include + `gamma_init()`)
- `src/info.c` (docstring)
- NEW `tests/test_gamma.c`; `tests/CMakeLists.txt` (COMMON_SRC + `gamma_tests`)
- `docs/spec/builtins/special-functions.md` (Gamma section)
- NEW `docs/spec/changelog/2026-06-08.md`; `Mathilda_spec.md` (changelog row)

Verification:
- `make -j` clean under `-std=c99 -Wall -Wextra`; no warnings on new files.
- `gamma_tests`: all 13 groups pass (exact ints/half-ints/BigInt, poles,
  infinities, symbolic, machine real/complex, MPFR precision tracking,
  incomplete + generalized, Listable, attributes).
- valgrind: a 16-evaluation Gamma driver matches the `Sin[1.0]` baseline
  byte-for-byte (12,800 B/400 blocks definitely + 3,720 B/56 indirectly —
  documented dyld/Accelerate noise); zero Gamma frames in any lost stack.
- REPL spot-checks reproduce every spec example value.

Design notes / deliberate limitations:
- Exact integer/half-integer `Gamma[z]` reuses the `Factorial[z-1]` machinery
  (DRY: one code path for exact int64 / BigInt / rational·Sqrt[Pi]).
- Real numerics use `tgamma` (machine) and `mpfr_gamma` / `mpfr_gamma_inc`
  (arbitrary). Machine complex uses a Lanczos approximation.
- Arbitrary-precision **complex** gamma is left symbolic on purpose — a fixed
  Lanczos series can't honour the advertised precision; honest > wrong.
- Out of scope (future): D/Integrate/Series of Gamma; closed forms for
  exact-integer incomplete gamma (`Gamma[2,3]` stays symbolic).
