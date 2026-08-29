# RischNormanBlake — task tracker

Engine: `src/calculus/int_rnb.c` (+ `int_rnb.h`), `Integrate`RischNormanBlake[f,x]`.
Scope: full algorithm incl. exact S-unit tier, n=1, K=Q(x). Reproduce §8 examples 1–4.

## Phase 0 — foundations / API study
- [ ] Read Expr construction/inspection API (expr.h), eval helpers, sym_names, internal_* calls
- [ ] Read intrischnorman.c helper region (mk_*, eval_* wrappers, driver, install idiom)
- [ ] Confirm how to call internal_factor / risch_squarefree_t / series / residue_compute / solve

## Phase 1 — scaffold + wiring (compile-green skeleton)
- [ ] int_rnb.h + int_rnb.c skeleton: builtin_rischnormanblake, int_rnb_init, install
- [ ] eval_* wrapper helpers (together/expand/numer/denom/cancel/gcd/coeff/degree/D)
- [ ] integrate.c wiring: include, init call, try_rischnormanblake, cascade insert, method enum/string/case, diagnostic+docstring lists
- [ ] Build clean (make -j), engine returns NULL for now (stub) — no regressions

## Phase 2 — RadicalField + element arithmetic
- [ ] radical detection in f -> (q, m); squarefree Qj; N1/N2 checks
- [ ] Ei, mul table, Lam; struct RadicalField
- [ ] rf_add/scal/mult/D/one/zero/mult_matrix/norm/inv/from_y_expr/to_y_expr/denominator
- [ ] micro-test element ops (mult table associativity, D(w_i)=Lam_i w_i)

## Phase 3 — heuristic parallel_integrate + exact degree bounds
- [ ] rnb_exact_degree_bounds (val_inf_w/val_inf_f, nu_i, e, g, Delta, k)
- [ ] rnb_parallel_integrate: ansatz, Dg, residual, coeff eqs, layered solve (mpq fast / symbolic)
- [ ] diff-back verify helper (flint_algebraic_field_normalize)
- [ ] Example 1 (1/Sqrt[x^2+1]) with cf-unit logand -> Log[x+y]

## Phase 3 — DONE: heuristic core works
- [x] f->element, denominator, ansatz, symbolic solve, diff-back verify
- [x] cases 1-4 (algebraic part only) verify to 0; non-elem/log cases decline cleanly

## Phase 4 — exact S-unit tier
- [x] cf_units (m=2 continued fraction) -> Example 1 = Log[x+Sqrt[1+x^2]] (check 0); c5 works
- [ ] places_over_denominator (ramified/unramified classification) + residues (debug surface)
- [ ] residue_at: ramified (trace) + unramified (branch series) — validate vs paper values
- [x] find_element_with_divisor (val conditions -> nullspace) + affine_divisor_ok (norm)
- [x] exact_logands driver (group residues -> divisors -> logands, torsion N search)
- [x] exact per-coordinate degree bounds (bounds.py) — tractable systems
- [x] KEY FIX: algebraic-constant abstraction (Sqrt->symbol) so Cancel/CoefficientList
      extract clean equations over Q(radicals); RootReduce eqs-satisfied verification
- [x] ALL FOUR paper examples verify (numerical diff-back = 0): 8.1 sqrt2/3, 8.2 torsion,
      8.3 omega, flagship dx/sqrt(x^2+1); + nested min case; 0.87s total

## Phase 5 — cascade, cleanup, tests, docs
- [ ] confirm Integrate[f,x] cascade routes radical integrands to RischNormanBlake
- [ ] clean output coefficients (RootReduce final gammas); silence Solve::svars/Power::infy
- [ ] remove RNB_DEBUG or keep getenv-guarded; fix misleading-indentation warning
- [ ] C unit test tests/test_int_rnb.c + stress corpus (.m) w/ numerical diff-back
- [ ] docstring done; docs/spec/builtins + changelog; make check-c99; regressions; valgrind

## Phase 5 — docs, tests, audits
- [ ] docstring; docs/spec/builtins + changelog
- [ ] tests/test_int_rnb.c (+ COMMON_SRC); .m REPL script
- [ ] make check-c99; regression suites; valgrind the 4 examples; graph refresh

## Review — COMPLETE

`Integrate`RischNormanBlake[f, x]` landed in `src/calculus/int_rnb.c` (+ `.h`),
wired into the Integrate cascade after `try_risch` and selectable via
`Method -> "RischNormanBlake"`. All four paper examples + a broad stress corpus
verify (numerical diff-back = 0); non-elementary/out-of-scope cases decline
cleanly (never wrong). C unit test `tests/test_int_rnb.c` passes; pmint,
integrals, and dispatch regressions unaffected; `make check-c99` clean.

Key implementation notes / gotchas encountered:
- **Style**: Expr*-orchestration over CAS wrappers, mirroring intrischnorman.c.
  Element = coordinate vector over the Trager basis w_i = y^i/E_i.
- **Leading-coefficient bug**: risch_squarefree_t returns MONIC Q_j, so the
  radicand's leading constant c = q/∏Q_j^j must ride the multiplication table on
  each m-th-power wraparound — without it `q = 1-x²` (c=-1) gave a sign-flipped
  (wrong) answer. Caught only by numerical diff-back, not by the algebraic
  eqs-satisfied check (which trusts the system it was handed).
- **Number-field linear algebra**: Mathilda's Cancel/Together/RootReduce do NOT
  cancel algebraic-constant common factors (Sqrt[3/2]/(Sqrt[3/2](x²-2)) stays;
  only the hang-prone FullSimplify reduces it). Fix: abstract each algebraic
  constant to a fresh SYMBOL before Cancel/CoefficientList (Cancel clears s/s),
  substitute back for Solve, verify equations with RootReduce.
- **Verification**: FullSimplify hangs on mixed poly+constant radical fields;
  use the extracted-equations RootReduce check + a high-precision numerical
  diff-back gate. The numerical gate is essential — it's the only check that
  catches a mis-built linear system (field-setup bug).
- **Robustness**: m≥3 residue tier is real-place only (complex-place branch
  series is a root-of-a-complex-number Puiseux expansion that hangs); plus a
  12s wall-clock budget. Pathological cases decline, never hang or mislead.
- **Leak**: fixed an unfreed `rows` buffer in rnb_find_element. Two residual
  one-time leaks are inside called library fns (risch_squarefree_t,
  builtin_linearsolve), not this engine.

Follow-ups (out of scope here): towers (n>1); complex-place m≥3 logands;
the per-call `Solve::svars` stderr line (Mathilda's Quiet is broken; the Cherry
engines emit it too); prettier (RootReduced) output coefficients.

## v0.122 — never hangs; complex-place integrands solve (2026-08-30)

Diagnosed the reported "struggles/slow/slow-to-fail" (all `Method ->
"RischNormanBlake"`): the hang was `Cancel`/`Together` on `I` + a symbolic
unknown in the parallel-solve accumulation (NOT the residue tier, which was
<0.2 s). Fixes in `src/calculus/int_rnb.c`:
- Linearise `Complex -> RNB$I0` before the solve, restore at equation extraction
  (sound: `I^2=-1`). Radicals are NOT abstracted pre-accumulation (`Sqrt[6]` vs
  `Sqrt[2]Sqrt[3]` would break the system) — only post-accumulation, as before.
  `(x^2-1)/((x^2+1)Sqrt[x^4+1])` now solves (conjugate complex logs).
- Hard wall-clock bounds: `TimeConstrained` caps on the heavy number-field calls
  (Series/NullSpace/Cancel/Together/Factor/LinearSolve/Coefficient) + deadline
  polls between the per-logand accumulation sub-ops; budget 4 s. No run hangs.
- `tests/test_int_rnb.c`: + `test_rnb_complex_places`; all groups pass.
- Acceptance (Method path, all 10 reported examples): 2 INTEGRATE (the working
  `Sqrt[x^4+1]` case + its complex-pole conjugate), 8 DECLINE-FAST (<5 s), 0 hang.

B3 follow-up (deferred deeper exact-tier rework — all deemed elementary, now
decline fast instead of hang): `(x-1)/((x+2)√(x^3-1))`,
`x/((x^3+8)√(x^3-1))`, `(x+1)/((x-1)√(x^4+x^2+1))`, `x/((2x^3-1)√(x^4-x))`,
`1/((x+1)(3x^2+1)^(1/3))`, `1/(x(x^2-3x+2)^(1/3))`,
`(1+2x^2)/((1+x^2)(1+3x^2)(2x^3-x)^(1/3))`, `1/(1+x^3)^(1/3)`. (Pre-existing,
non-RNB: the default cascade still hangs on some of these via a later engine.)
