# Reduce — implementation tracker

Design: `REDUCE_PLAN.md` (repo root). Full CAD roadmap + full companion family.

## Phase 0 — Front-end + logical normal-form skeleton  (DONE 2026-08-23)

- [x] `src/solve/reduce_form.{h,c}` — RRel / RAtom / RConj / RForm (DNF), builders,
      `rform_or`/`rform_and`/`ratom_negate`, `rform_simplify`, `rform_to_expr`.
- [x] `src/solve/reduce_atom.c` — `reduce_atom_canonicalize` (poly REL 0), constant-atom
      decide via `evaluate`, `reduce_form_from_expr` (parse &&/||/!/Implies/Xor/Inequality → DNF).
- [x] `src/solve/reduce.{h,c}` — `builtin_reduce` (Solve-style parse, True/False,
      bad-var warn, dom positional), dispatch skeleton (constant or NULL), `reduce_init`.
- [x] `src/sym_names.{h,c}` — `SYM_Reduce` (3 sites).
- [x] `src/core.c` — `#include "reduce.h"` + `reduce_init()`.
- [x] `tests/test_reduce.c` + `tests/CMakeLists.txt` (COMMON_SRC + `reduce_tests` + `add_test`).
- [x] docs + changelog.
- [x] main build clean (gcc-16); `reduce_tests` 24/24 pass; `solve_tests` still green;
      `make check-c99` clean; `leaks` 0 leaks.

Verified: `Reduce[True,x]→True`; `Reduce[1<2,x]→True`; `Reduce[x==x,x]→True`;
`Reduce[3<2,x]→False`; `Reduce[x>0,x]`/`Reduce[x^2==4,x]` stay unevaluated;
`Reduce[x==1,5]` → `Reduce::ivar` + unevaluated.

## Phase 1 — Complete univariate equation solver (Complexes)  (DONE 2026-08-23)

- [x] `src/solve/reduce_eq.{c,h}` — `reduce_eq_univariate`, lc-vanishing recursion;
      generic roots via `Solve`; degree/coeff via `Exponent`/`Coefficient`/`Expand`.
- [x] DNF-layer additions: `RAtom.display` (solved-form emission), `ratom_solved`,
      sign-normalization of EQ/NE atoms (`-b==0` → `b==0`).
- [x] `reduce.c` routing: single univariate EQ atom over Complexes → equation engine.
- [x] CMake COMMON_SRC += reduce_eq.c; tests extended (`test_equations`).
- [x] docs + changelog.
- [x] `reduce_tests` 32/32 pass; `solve_tests` green; `check-c99` clean; `leaks` 0.

Verified: `a x==b → (a!=0 && x==b/a) || (a==0 && b==0)`; `x^2==4 → x==-2||x==2`;
`x^2==-1 → x==-I||x==I`; `2x==6 → x==3`; `a x^2+b x+c==0` → full 3-level split.

## Phase 2 — Univariate real sign diagram (Reals)  (DONE 2026-08-23)

- [x] `src/solve/reduce_univar.{c,h}` — sign diagram: roots via `Solve[..,Reals]`,
      breakpoint sort/dedup + interval sampling, exact sign via native-rational +
      `flint_qqbar_compare`, union-of-cells emission (Inequality chains, cofinite `!=`).
- [x] `reduce.c` routing: `reals && nv==1` → `reduce_univar` (any poly eq/ineq combo).
- [x] CMake COMMON_SRC += reduce_univar.c; tests extended (`test_real_inequalities`).
- [x] docs + changelog.
- [x] `reduce_tests` all pass (41 assertions); `solve_tests` green; `check-c99` clean;
      `leaks` 0; main build warning-free.

Verified: `x^2>1→x<-1||x>1`; `x^2>=1→x<=-1||x>=1`; `x^2<1→-1<x<1`;
`(x-1)(x-2)(x-3)>0→1<x<2||x>3`; `x^2!=1→x!=-1&&x!=1`; `x^2==4→x==-2||x==2`;
`x^2<2→-Sqrt[2]<x<Sqrt[2]`; `x^2+1>0→True`; `x^2+1<0→False`; `x>0&&x<1→0<x<1`.

## Phase 3 — Linear real systems via Fourier–Motzkin (Reals)  (DONE 2026-08-23)

- [x] `src/solve/reduce_fm.{c,h}` — FM elimination over exact-rational coeff vectors:
      project vars last→first, feasibility from the fully-projected constants,
      triangular emission (bounds per var, `==` re-detection, free-var omission),
      DNF handled conjunct-by-conjunct then OR-ed.
- [x] `reduce.c` routing: `reals && nv>=2` → `reduce_fm` (declines if non-linear).
- [x] **Soundness fix**: `RAtom.nonconst_denom` flag — inequalities that cleared a
      variable denominator (`1/x<1`) are no longer constant-decided or handled by the
      real engines; Reduce declines instead of answering wrong. (Fixed a Phase-2 bug.)
- [x] CMake COMMON_SRC += reduce_fm.c; tests extended (`test_linear_systems`).
- [x] docs + changelog (Phase 3 + soundness fix).
- [x] `reduce_tests` all pass; `solve_tests` green; `check-c99` clean; `leaks` 0;
      build warning-free.

Verified: `x+y<1&&x>0&&y>0 → 0<x<1 && 0<y<1-x`; `x+y==1&&x>0 → x>0 && y==1-x`;
`x>1&&x<0&&y>0 → False`; `2x+3y<=6&&x>=0&&y>=0 → 0<=x<=3 && 0<=y<=2-2x/3`;
`x<0||x>1 → x<0||x>1`; nonlinear/`1/x<1` decline.

## Phase 5 — Integers / Rationals (thin wrapper)  (DONE 2026-08-23)

- [x] `src/solve/reduce_int.{c,h}` — reformat `Solve[..,dom]` output into logical
      form (`||` of `==` atoms; `Element[C[k],dom]` for parametric families).
- [x] `reduce_univar_integers` (in reduce_univar.c) — bounded integer enumeration
      over the sign diagram for the inequality case Solve declines; shared
      `collect_breakpoints` factored out of reduce_univar.
- [x] `reduce.c` routing: `dom ∈ {Integers, Rationals}` → `reduce_integers`.
- [x] CMake COMMON_SRC += reduce_int.c; tests extended (`test_integer_domain`).
- [x] docs + changelog.
- [x] `reduce_tests` all pass; `solve_tests` green; `check-c99` clean; `leaks` 0;
      build warning-free.

Verified: `x^2==4→x==-2||x==2`; `x^2<10&&x>0→x==1||x==2||x==3`; `1<=x<=3→...`;
`x+y==5&&x>0&&y>0→4 tuples`; `2x+3y==1→C[1]∈Integers && x==-1+3C[1] && y==1-2C[1]`;
`x^2==2→False`; Rationals; unbounded `x>0` declines.

## Phase 4 — Parametric linear systems over Complexes  (DONE 2026-08-23)

- [x] `src/solve/reduce_sys.{c,h}` — symbolic Gaussian elimination with case
      splitting (nonzero-const pivot direct; symbolic pivot p → p!=0 branch +
      p==0 branch via Solve-substitute-recurse); back-substitution (graft) for
      per-variable param expressions; LSol DNF-of-cases intermediate.
- [x] `reduce.c` routing: complexes && all-EQ → Phase 1 (single univar eq) else
      Phase 4 (`reduce_eq_system`, declines if non-linear).
- [x] **Bug fixed**: double-free of `prod` in the elimination loop (it is consumed
      by the `Subtract` node) — corrupted the heap → runaway evaluator recursion.
- [x] CMake COMMON_SRC += reduce_sys.c; tests extended (`test_parametric_systems`,
      test_unevaluated's stale `x+y==1` case replaced with a nonlinear system).
- [x] `reduce_tests` all pass; `solve_tests` green; `check-c99` clean; `leaks` 0;
      valgrind: no code-level memory errors; build warning-free.

Verified: `a x+y==1 && x+y==0 → 1-a!=0 && x==1/(a-1) && y==1/(1-a)`;
`a x==1 && x==2 → 2a-1==0 && x==2`; `x+y==1 → x==1-y` (free var);
`x+y==3 && x-y==1 → y==1 && x==2`; 3-var system; nonlinear system declines.

Next: **Phase 6** (multivariate nonlinear CAD — the big one), **7** (QE:
Exists/ForAll/Resolve), **8** (companions: LogicalExpand/FindInstance/
CylindricalDecomposition + polish: default-domain inequality → Reals,
unbounded-integer inequality forms, nicer parametric-condition display).

## Later phases
1 Complete univariate equations · 2 Univariate real sign diagram · 3 Fourier–Motzkin ·
4 Parametric linear systems · 5 Integers/Rationals · 6 CAD · 7 QE · 8 Companions + polish.
