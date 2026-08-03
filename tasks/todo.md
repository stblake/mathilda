# Task: Boolean NDArray dtype + Compile support

Full plan: `~/.claude/plans/let-s-extend-ndarray-to-memoized-stonebraker.md`

## Plan (checklist)

- [x] **Phase 1** — `NDT_BOOL` storage dtype (`src/expr.h`, `src/ndarray.c`): the
  ~13 `ndt_*` choke points, string↔enum, element→Expr → `True`/`False`,
  `leaf_to_component`, constructor error/docstring. Numeric engines DECLINE bool
  (delist to symbolic).
- [x] **Phase 2** — auto-packing of `True`/`False` lists (`src/pack.c`: `PK_BOOL`).
- [x] **Phase 3** — producers/consumers (gap C.1): sign predicates → packed bool
  array (`ndint_sign_predicate`), `Boole` → int64, `AllTrue`/`AnyTrue`/`NoneTrue`
  byte-scan (`funcprog.c`). `Boole` added to `pack.c` AWARE.
- [x] **Phase 4** — `Compile[]` bool arrays: `SYM_Boolean`, argspec, boundary
  marshalling, `A_LOAD_B`/`A_STORE_B` opcodes, `ct_elem_ndt`, `ct_is_elem`.
- [x] **Phase 5** — tests, docs, audit baselines.

## Review

**What shipped.** A `"bool"` NDArray dtype (`"Boolean"` alias), 1 byte/element,
storable/printable/constructible/auto-packable. Sign predicates now emit packed
bool arrays (closing performance.md §13 gap C.1); `Boole`/`AllTrue`/`AnyTrue`/
`NoneTrue` consume them off the buffer. `Compile[]` takes/returns/indexes bool
arrays and declares `_Boolean` scalars (which were internally `CT_BOOL` but had
no argspec spelling until now).

**Two load-bearing design decisions.**
1. *Bool is not numeric.* Arithmetic/transcendentals over a bool array return
   `NULL` from every numeric engine and delist to the symbolic `List`
   (`Sin[True]`, `True + True` unevaluated — Mathematica-faithful; mirrors the
   compiler's existing "Bool is not numeric" `num_common` invariant).
   `ndarray_int64_delist_retry` was broadened to also delist bool so Plus/Times/
   Power thread correctly.
2. *Bool producer results present as a `List`, not the input's presentation.* The
   sign predicates have always answered with a `List`; `test_ndarray_functions.c`
   pins it. So the result is a PACKED bool List (observably identical, now on a
   buffer), unlike `Sin` which inherits its input's presentation. Same for `Boole`.

**Verification.** Clean build (`-std=c99 -O3`), `make check-c99` green.
`check-compile-coverage` green (sign predicates moved EXEMPT→BASELINE — a bool
lowering is now *possible*, only undone). `check-packed-aware` green
(`Boole` added to AWARE). `check-array-exactness` green on the final binary
(0 MIXED). 21+ test suites pass, including new bool cases in `test_ndarray.c`,
`test_packed_list.c`, `test_compiledfunction.c`. `check-nd-surfaces` /
`check-fastpath-sweep` deliberately NOT run (slow; per project guidance).

**Deferred / non-goals (as agreed).** Comparison operators do not thread over
arrays (not Listable in WL). No fused whole-array logical ops in Compile. No bool
arithmetic. A Compile *array* lowering for the sign predicates is now possible but
undone (tracked in `compile_coverage.py` BASELINE).
