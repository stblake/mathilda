# Task: Make Interpolation / ListInterpolation build as efficient as possible

Plan: `/Users/user/.claude/plans/peaceful-kindling-adleman.md`

Root cause (profiled): `InterpolatingFunction` stores `table` as a boxed List of
n `{coord,value}` pairs → build churns ~3n Expr nodes (~28–38 ms/10⁵). The numeric
build itself is ~2 ms (NO_PACK). Fix: store the table as a packed n×(m+1) float64
NDArray for machine-real scalar-valued interpolants.

## Phase 0 — infra / recon
- [ ] Read NDArray struct + constructor + present_as (transparent packed List)
- [ ] Read IFun struct (add a `scalar_valued`/`Vfilled` flag + own V buffer)
- [ ] Confirm `is_ndarray` predicate + buffer layout (row-major)

## Phase 1 — packed-table helpers (src/interp.c)
- [ ] `interp_table_is_packed(table)` predicate
- [ ] Packed accessors: coord(i,k), value(i), npts, m from a packed table
- [ ] Builder: doubles[] (n×(m+1)) → transparent packed NDArray Expr

## Phase 2 — build sites (the win)
- [ ] `builtin_interpolation_impl`: gate (scalar-valued ∧ Ksupplied==0 ∧ machine-real) + build packed table
- [ ] `builtin_listinterpolation`: build packed table straight from value buffer + xs grids
- [ ] `builtin_interpolation`: route real NDArray input into packed build (no delist)

## Phase 3 — readers (packed branch each)
- [ ] `build_grid`: coords per column; fill `f->V` directly; skip `entryAt`
- [ ] `table_Ksupplied` → 0; `obj_cache_load` → machine, non-mpfr
- [ ] exact-node query (interp_apply 1464–1481)
- [ ] `interp_eval_double` / `interp_vector_1d`: honour scalar_valued flag
- [ ] `integrate_interp.c`: packed branch or delist-fallback

## Phase 4 — verify
- [ ] `tests/test_interp.c`: packed≡boxed equivalence (1-D/2-D/Spline/Hermite/order1/vector) + exact-node + exact/MPFR/array-valued stay boxed
- [ ] `make` clean; `make check-c99`
- [ ] surface audits: check-nd-surfaces / check-array-exactness / check-packed-aware / check-compile-coverage
- [ ] valgrind leak-clean on 10⁵ build+eval
- [ ] benchmark 16 re-run: build rows ≤1.5×, CHECK-FAIL=0, eval unchanged
- [ ] changelog + docs

## Review

**Outcome: all four `Interpolation`/`ListInterpolation` build rows moved from
4.5–13× behind scipy to AHEAD.** (ListInterp 1-D 30→0.82ms; 2-D 17→0.34ms;
Interp[NDArray] 34→0.42ms.) Eval rows unchanged (still parity/ahead).

Implementation (all in `src/interp.c` unless noted):
- Packed `n×(m+1)` float64 table representation for machine-real, scalar-valued,
  value-only (`Ksupplied==0`, non-Hermite) interpolants. Stored as a VISIBLE
  NDArray (not a packed List — a packed List would be materialised by the
  transparency gate into flat `{c,v}` triples the boxed readers can't parse).
- Producers: `interp_try_pack_data` (impl gate), `builtin_listinterpolation`
  direct packed build with double-only abscissae (`listinterp_axis_double` +
  `listinterp_pack_emit`/`_buffer`), and `interp_ndarray_fast` (buffer-direct for
  `Interpolation[NDArray]`, no delist). Shared `interp_finish_object` /
  `interp_parse_option` / `interp_resolve_periodic` helpers (impl refactored onto
  them too).
- Readers: `build_grid_packed` (fills `f->V` directly, `entryAt=NULL`);
  packed branches in `table_Ksupplied`, `obj_cache_load`, the exact-node query,
  and the two `value_shape` sites.
- `src/calculus/integrate_interp.c`: delist-fallback so `Integrate[ifun]` (a cold
  path) still reads a packed table (regression fixed + covered).

Verification: interp/interp_poly/ndsolve suites green; `test_packed_table` added
(8 assertions incl. Integrate + exact-Integer-stays-boxed). check-c99,
check-packed-aware, check-array-exactness (0 MIXED), check-nd-surfaces all green.
valgrind: definitely-lost identical at 1× and 12× (13,440B baseline — no per-call
leak). check-compile-coverage: exempted the 3 interpolation heads (they return
InterpolatingFunction/polynomial objects); the gate has a large PRE-EXISTING
unrelated backlog (ArrayPlot/Fit/GaussianFilter/… — already red on main, not
touched here).
