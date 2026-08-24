# Fix: bare `N[expr]` must target machine precision

## Bug
`N[N[Pi, 100]]` returned a 100-digit number (precision 100.243); Mathematica
returns a machine-precision number. `N[N[Pi, 100], 30]` already worked (30.103).
Root cause: the one-argument `N` builtin set `spec.preserve_inexact = true`,
which told `numericalize` to keep an already-approximate MPFR leaf at its
existing precision under MACHINE mode.

## Decision
The `preserve_inexact` field exists ONLY to implement this wrong behavior:
every call site except `builtin_n` sets it `false`, so its three consult sites
always took the false path everywhere else. Remove the field at the root.

## Steps
- [ ] `src/numeric.c` EXPR_MPFR/MACHINE branch — always down-convert (drop flag branch)
- [ ] `src/numeric.c` `builtin_n` one-arg — stop setting the flag; fix comment
- [ ] `src/numeric.c` `numeric_plan_working_spec` — drop the flag bail + `= false`
- [ ] `src/numeric.c` top-level `numericalize` inf/zero guard — drop flag term + `= false`
- [ ] `src/numeric.c` `parse_precision_arg` — drop `= false`
- [ ] `src/numeric.h` — remove field, `numeric_machine_spec` init, fix stale comment
- [ ] External call sites (nsum, nlimit, nderiv, nresidue, nroots, nseries,
      ndsolve, piecewise, random, root_numeric) — drop `= false` sets
- [ ] Book: rewrite note `arithmetic/contagion/2`; regenerate transcript
- [ ] Build + verify N[N[Pi,100]] machine, N[N[Pi,100],30]=30, N[N[Pi,100],200] stays 100
- [ ] Run numeric tests; docs/changelog; memory update

## Review — DONE

Root-caused to `builtin_n` setting `NumericSpec.preserve_inexact = true` for the
one-argument form. That flag existed only to implement the wrong behavior — every
other of its ~18 sites set it `false`, so its three consult branches always took
the false path elsewhere. Removed the field entirely (behavioral no-op except at
`builtin_n`), plus the now-dead `ExactScan.has_mpfr_leaf` it was the only reader of.

Verified (fixed binary):
- `N[N[Pi,100]]` → `3.14159`, `Precision` → `MachinePrecision` (was 100.243)
- `N[N[Pi,100],30]` → 30.103 (unchanged); `N[N[Pi,100],200]` → stays 100.243
- Contagion `1.+N[Pi,100]` → machine; `SetPrecision[Pi,50]` → 50.272; `N[E,50]` →
  50.272; `N[Exp[1000]]` → finite `1.97e+434` (inf/zero MPFR fallback intact)

Tests: rewrote `test_numeric.c::test_n_bare_targets_machine_precision` and the
nested-N asserts in `test_numeric_largearg.c`; both suites + nsum/nlimit/nderiv/
nresidue/nseries/nroots/root_numeric/piecewise/accuracygoal/machine_number_q all
pass. `make check-c99` clean. Clean `-Wall -Wextra` build.

Docs/book/memory: `docs/spec/builtins/arithmetic.md`, weekly changelog,
`book/chapters/math/arithmetic.tex` note + regenerated transcript
(`Out[2]= 3.14159`), memory `project_n_preserve_inexact` rewritten + MEMORY.md
hook, lesson appended to `tasks/lessons.md`.
