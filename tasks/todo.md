# Issue #52 — iterator termination in double precision (+ Table cap)

## Root cause
Table/Do/Sum/Product drive loop termination on a `double val <= double max_val`
test. Near 2^63, consecutive int64 values collapse to the same `double` (ULP of
2^63 is 2048 > a unit step), so the test never fires:
- **Do** infinite-loops; **Sum/Product** run to their 10^8 term cap; **Table**
  runs to a **10^6** cap and *silently truncates*.
The exact running value `curr_e` is maintained alongside but never consulted for
termination. Separately, Table's 10^6 cap truncates *legitimate* tables > 1M
(`Table[i,{i,1,2000000}]` -> 1000001), and truncates silently unlike Sum/Product
which return unevaluated.

## Fix
- [x] `iter_range_continue()` shared helper (iter.c/iter.h): int64 fast path
      (exact, no GMP — preserves PR #50 loop perf), GMP compare when a BigInt is
      involved, double test for real / rational-exact, `is_inf` short-circuit.
- [x] Do — replace while-condition (iter.c).
- [x] Table — replace fill + pre-count conditions; raise cap 10^6 -> 10^8;
      return unevaluated on overflow (align with Sum/Product).
- [x] Sum — replace loop condition (sum.c expand_range).
- [x] Product — replace loop condition (product.c expand_range).
- [x] Auto-compile: numloop_do_range's three int64 loops advance with the
      overflow-checked ci_add_i64 (stop at the edge instead of wrapping).
- [x] Compile[]: default mode already bails to the (now-fixed) interpreter;
      wrap mode ("Speed") kept loop-control arithmetic checked via new
      IF_FORCECHK bit (compile_internal.h + compile.c funnel + 3 emit sites).
- [x] Perf: iter_range_continue inlined (BigInt arm out-of-line) so the tight
      loop is flat vs baseline (0.0081 vs 0.0079 on Do[Null,{i,200000}]).
- [x] Extensive unit tests: 6 fns in test_iter.c + test_cf_loop_counter_boundary.
- [x] Build clean, check-c99 clean, 13-suite regression sweep green, changelog + docs.

## Review
Fixed across all three layers (interpreter, auto-compile, Compile[]) plus the
separate Table 10^6-truncation bug the investigation surfaced. Verified:
- `Table[i,{i,9223372036854775805,9223372036854775807}]` -> 3 elements ✓
- `Do`/`Sum`/`Product` same span terminate ✓
- `Table[i,{i,1,2000000}]` -> 2,000,000 (cap no longer truncates) ✓
- over-cap exact range -> unevaluated instantly ✓
- Compile default + "Speed" boundary -> correct (was hang / under-run) ✓
- exactness unchanged; no perf regression ✓
