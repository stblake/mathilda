# RankedMin / RankedMax — the n-th order statistic

Plan: `~/.claude/plans/let-s-implement-rankedmin-and-silly-umbrella.md`

Core identity: `RankedMax[list, k] ≡ RankedMin[list, -k]`. Select ascending rank `r`.

## Tasks
- [x] 1. Symbols in `sym_names.c` (defs + init) and `sym_names.h` (extern)
- [x] 2. NDArray buffer fast path `ndred_ranked_min/max` in `ndreduce.c` + `ndreduce.h`
- [x] 3. Builtins `builtin_ranked_min/max` + `ranked_select` (numeric-key quickselect) in `sort.c` + `sort.h`
- [x] 4. Register in `core.c` (Protected)
- [x] 5. Packing: AWARE + INT64_OK in `pack.c`
- [x] 6. Compile lowering: `compile.c` (NdRedSpec.nextra, lookup, emit, infer, VM) + `compile_internal.h` (V_NDREDN) + `disasm.c`
- [x] 7. Docstrings in `info.c`
- [x] 8. Tests: new `test_ranked.c` (9 groups) + CMake; 5 compile parity cases in `test_compile.c`
- [x] 9. Docs: `structural-manipulation.md` + changelog `2026-08-03.md`
- [x] 10. Verify: build, tests, check-compile-coverage, check-c99, valgrind

## Review

**What shipped.** `RankedMin[list, n]` / `RankedMax[list, n]` — the n-th smallest /
largest element (negative n counts from the other end). One routine: `RankedMax`
negates its rank and both reduce to selecting an ascending rank `r`. So
`RankedMin[l, 1]`==`Min[l]`, `RankedMin[l, -1]`==`Max[l]`. Interpreter + packed
buffer fast path + `Compile[]` + auto-compilation.

**Design — two ordering paths (sort.c).** A definite result needs every element to
be a real number. One scan chooses: an all-`is_real_numeric` list orders by
`expr_compare` (BigInt/Rational-exact); otherwise a machine-double key orders
symbolic reals (`Pi`, `E`, `Sqrt[2]`, `Pi+E` via `numericalize`) and `±Infinity`
(→`±HUGE_VAL`), bailing to NULL on a free symbol / non-real complex. Both paths
return the element in its **exact** form and select via an O(m) quickselect over an
index array (comparator = key/expr_compare, original index stable tiebreak,
`cmp(i,i)==0`). High-precision reals stay untouched — only their double key is taken.

**Compile — first `(array, int) → scalar` lowering.** No existing table fit
(`ND_REDS` was `na==1`-gated, `ND_FNS` returns an array). Added `NdRedSpec.nextra`
(existing rows zero-fill to 0), a new opcode `V_NDREDN` = `V_NDRED`'s scalar write +
`A_NDFN`'s trailing-int read (c->a=array, c->b=n), and `nd_red_lookup` now keys on
`na == 1 + nextra`. Emitted via `arr_op(c, OP_V_NDREDN, a_array, n_val, ...)` which
wires both registers and recycles both temps. `CompileDiagnostics[…, RankedMin[v,k]]`
→ `Compiled -> True` at real & int shapes.

**Buffer path (ndreduce.c).** `ndred_ranked_min/max` mirror `ndred_median`: int64
selects exactly via `nd_sort_i64_asc` (Integer, not a rounded Real), real via
`nd_select_kth` (O(m)); complex / rank>1 / out-of-range degrade to
`ndarray_delist_and_reeval`. Shared by the interpreter builtin and the Compile VM.

**Verification.** Clean `-O3 -std=c99` build. `test_ranked.c`: all 9 groups pass
(every documented example, the identities, infinities, exact BigInt/Rational,
duplicates, symbolic reals, unevaluated cases, NDArray int64/float64). 5 compile
parity cases exact (max_rel=0.0). `check-c99` green; `check-compile-coverage` green
(no new gaps — RankedMin/Max compile at `v_rk`/`v_ik`). sort/stats/compile/coverage
suites pass. Valgrind: leak profile **byte-identical** to a matched Min/Max/Sort
baseline (13,376 B / 418 blocks — the macOS objc baseline), none of my functions in
any leak stack → zero new leaks.
