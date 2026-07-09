# Rename Matrix[...] → NDArray[...] + N-dimensional (numpy conventions)

Decisions: full internal rename; introspection now, broadcasting/higher-rank Dot later.

## Phase 1 — Full rename (Matrix → NDArray)
- [ ] git mv src/matrix.c → src/ndarray.c, src/matrix.h → src/ndarray.h
- [ ] Rewrite ndarray.h/.c: guard MATHILDA_NDARRAY_H, fns ndarray_* , builtin_ndarray,
      ndarray_init, register "NDArray", NDARRAY_MAX_RANK, docstring
- [ ] expr.h: EXPR_MATRIX→EXPR_NDARRAY, MatrixData→NDArrayData, .matrix→.ndarray,
      expr_new_matrix→expr_new_ndarray
- [ ] expr.c: same identifiers
- [ ] Targeted edits (EXPR_MATRIX / data.matrix / API fns / include): eval.c, match.c,
      numeric.c, precision.c, print.c, part.c, sort.c, simp/simp_complexity.c, core.c,
      plus.c, times.c, linalg/dot.c, list/listpredicates.c, calculus/series.c
- [ ] sym_names.{c,h}: SYM_Matrix→SYM_NDArray, intern "NDArray"
- [ ] Do NOT touch: MatrixQ / MatrixPower / MatrixRank / MatrixForm / *Matrix builtins,
      stats.c `is_matrix` local, vm_is_matrix

## Phase 2 — N-dim correctness + numpy introspection
- [ ] Fix Length[NDArray] → shape[0] (numpy len; currently 0)
- [ ] Add ArrayDepth (ndim for NDArray; rectangular depth for lists) — reuse get_dimensions
- [ ] Add ArrayQ (True for NDArray; rectangular full-depth list)
- [ ] Keep VectorQ(rank 1)/MatrixQ(rank 2); Dimensions=shape (exists)

## Phase 3 — docs / tests / changelog
- [ ] test_matrix.c → test_ndarray.c; rename symbols; add rank-3/4 cases
- [ ] tests/CMakeLists.txt: src/matrix.c→src/ndarray.c; matrix_tests→ndarray_tests (3 sites)
- [ ] docs/spec/builtins/linear-algebra.md: ## Matrix → ## NDArray (ndim, ArrayDepth, ArrayQ)
- [ ] docs/spec/changelog/2026-07-06.md: rename + ndim entry
- [ ] docstrings + attributes for NDArray/ArrayDepth/ArrayQ

## Phase 4 — build + verify
- [ ] make clean && make (clean: stale-object ABI trap on expr.h change)
- [ ] ndarray_tests pass
- [ ] Smoke: rank 1/2/3/4 construct, Dimensions, ArrayDepth, ArrayQ, Length, Head,
      MatrixQ/VectorQ, Dot, Plus/Times elementwise, Normal, Precision/Accuracy/MatchQ

## Deferred (follow-up, noted in changelog)
- numpy broadcasting (scalar+array, shape-compatible elementwise)
- higher-rank Dot / matmul batching

## Review (done 2026-07-09)
- Full internal rename Matrix→NDArray: EXPR_NDARRAY, NDArrayData, data.ndarray,
  src/ndarray.{c,h}, ndarray_*, SYM_NDArray. MatrixQ/MatrixPower/MatrixRank etc.
  untouched. Two printer literals ("Matrix[" fullform, "\text{Matrix}" tex) were
  string-literals the token-sed missed — fixed by hand (user caught one).
- Per user course-corrections: dropped ArrayDepth (use Depth[NDArray]=rank+1);
  ArrayQ→NDArrayQ (simple type predicate); no Part indexing (over-design);
  Length[NDArray]=shape[0] kept.
- NDArray::shape warning on incompatible-shape Plus/Times (all-NDArray operands),
  leave unevaluated like Dot::dotsh; deduped per-eval via eval clock.
- Verified: clean build (0 warnings), 29 ndarray_tests, valgrind clean (no
  NDArray frames), full functional battery. assert_eval_eq is NDEBUG-hollow —
  real verification via pipe protocol.
