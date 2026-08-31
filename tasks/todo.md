# DSolve — Riccati (1a)

`y'[x] == q0(x) + q1(x) y + q2(x) y^2` (q2 != 0), linearised by `y = -u'/(q2 u)`
to the 2nd-order linear ODE `u'' - (q1 + q2'/q2) u' + (q0 q2) u == 0`, solved by
recursing into the scalar cascade (M5 engine), then mapped back. One constant
(collapse `C[2] -> 1`).

## Implementation
- [ ] New `src/calculus/dsolve_riccati.c` (try + builtin + init), mirroring
      dsolve_bernoulli.c (coeff extraction) + dsolve_reduce_order.c (recursion)
- [ ] Wire into `src/calculus/dsolve.c`: enum `DS_RICCATI`, string map, externs,
      cascade line after `dsolve_fos_try`, pinned case, `dsolve_riccati_init()`
- [ ] Add `../src/calculus/dsolve_riccati.c` to `mathilda_common` in tests/CMakeLists.txt

## Tests
- [ ] tests/test_dsolve.c: t_method_riccati, t_riccati_more (incl. y'==y^2+x Airy),
      t_ivp_riccati; register in main()
- [ ] tests/test_dsolve_stress.c: riccati_ok forward generator (from-spectrum) +
      t_stress_riccati; register in main()

## Docs
- [ ] docs/spec/builtins/calculus.md: DSolve`Riccati methods-table row
- [ ] docs/spec/changelog/2026-08-31.md: "DSolve — Riccati (1a)" section
- [ ] DSOLVE_PLAN.md §1a: flip [ ] Riccati -> [✓]

## Gates
- [ ] make -j (main binary) clean
- [ ] REPL spot-checks (y'==y^2+x; from-spectrum; IVP; fos unchanged; ?DSolve`Riccati)
- [ ] ctest -R dsolve (3/3 green)
- [ ] make check-c99 (exit 0)
- [ ] valgrind: one-time engine baseline, no per-call delta
- [ ] rebuild code-review graph

## Review

Done. Riccati (§1a) implemented as one new file + wiring + docs + tests, reusing
the M5 2nd-order linear engine via the classical `y = -u'/(q2 u)` linearisation.

- `src/calculus/dsolve_riccati.c` — coeff extraction (Bernoulli pattern) +
  recurse-into-scalar-engine (ReductionOfOrder pattern); `C[2]->1` collapse to the
  single Riccati parameter. Cascade slot after `fos`, before `autonomous`.
- Solves the Airy-linearised `y'==y^2+x` (previously declined), the elementary
  from-spectrum family, variable-coefficient `y'==x+x y^2`, and IVPs; the pinned
  method declines on a linear equation (`q2==0`). `fos` still owns `(x+y)^2`.
- Tests: 3 unit (`t_method_riccati`, `t_riccati_more`, `t_ivp_riccati`) + a
  from-spectrum stress family (`t_stress_riccati`, 7 spectra).

Gates: `ctest -R dsolve` 3/3 green (18.0 + 6.1 + 8.5 s); `make check-c99` exit 0;
valgrind 1× == 6× = 13,440 def + 6,312 indir (one-time engine baseline, no
per-call leak); code-review graph rebuilt.

No user correction occurred; no lesson added. LSP (clangd) noise about
`ATTR_READPROTECTED` / missing includes is stale — the GCC/CMake build is clean.
