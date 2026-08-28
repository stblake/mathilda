# Todo: Reduce Phase 8 — CylindricalDecomposition

Plan: `/Users/user/.claude/plans/let-s-continue-the-implementation-robust-acorn.md`

## Implementation
- [ ] `SYM_CylindricalDecomposition` at 3 sites (`sym_names.h`, `sym_names.c` decl+intern)
- [ ] `builtin_cylindrical_decomposition` in `reduce_companions.c` (Reals-forcing delegate to Reduce)
- [ ] Register in `reduce_companions_init`: builtin + ATTR_PROTECTED + docstring
- [ ] Prototype + header-comment update in `reduce_companions.h`
- [ ] Bump `src/version.h` 0.110 → 0.111

## Tests / verify
- [ ] `test_cylindrical_decomposition` in `tests/test_reduce.c` + TEST() registration
- [ ] Build; probe behaviour via NDJSON pipe; capture FullForm
- [ ] Run reduce_tests + reduce corpus (expect 160/160)
- [ ] `make check-c99` clean; valgrind at macOS baseline
- [ ] Audits: check-packed-aware / nd-surfaces / fastpath-sweep / compile-coverage (EXEMPT only if flagged)

## Docs
- [ ] `docs/spec/builtins/solutions-of-equations.md` — new `## CylindricalDecomposition`
- [ ] `docs/spec/changelog/2026-08-24.md` — dated section
- [ ] `REDUCE_PLAN.md` — Phase-8 status flip

## Review — complete & verified

**`CylindricalDecomposition` (Phase 8, v0.111) shipped.** A Reals-only front-end
(`src/solve/reduce_companions.c`, `builtin_cylindrical_decomposition`) that validates
arity (`[expr,vars]`, or redundant `[expr,vars,Reals]`; other 3rd positional declines) and
`vars` (symbol or List of symbols), then builds+evaluates `Reduce[expr, vars, Reals,
<trailing option Rules…>]` under message suppression and declines iff the result is still
headed by `Reduce`. Zero engine duplication — the whole `Reduce` pipeline (preprocessing,
DNF build, per-arity Reals dispatch: FM / sign diagram / CAD) is reused. Merged cylindrical
output, Mathematica-faithful.

- **Registration:** `SYM_CylindricalDecomposition` (3 sites), `symtab_add_builtin` +
  `ATTR_PROTECTED` + docstring in `reduce_companions_init`. Version 0.110 → 0.111.
- **Tests:** `test_cylindrical_decomposition` (14 asserts) — cylindrical forms, True/False,
  list-as-conjunction, bare-symbol var, and 3 sound declines. `reduce_tests` all pass;
  reduce corpus 160/160; Reduce/FindInstance unchanged (regression probe).
- **Portability/memory:** `make check-c99` clean; valgrind 13,440/420 = macOS baseline,
  **no new leak**. `check-packed-aware` OK (CD not flagged — symbolic head).
- **Docs:** spec section `## CylindricalDecomposition`; changelog `2026-08-24.md`;
  `REDUCE_PLAN.md` Phase-8 status flipped to ✅ + dated deviation note.

**Flagged (pre-existing, out of scope):** `make check-compile-coverage` fails on 22 NEW
heads (`Image*`, `Interpolation`, `Fit`, `GaussianFilter`, …) that gained numeric fast paths
without a BASELINE/lowering entry — unrelated to this change (none in the diff; CD is a
symbolic head and is not flagged). Belongs to whoever added those fast paths.
