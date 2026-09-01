# ListGradient — implementation todo

## Core module
- [x] `src/list/list_gradient.h` — header + `list_gradient_ndarray` proto
- [x] `src/list/list_gradient.c`
  - [x] Fornberg weights: `fd_weights_double`, `fd_weights_expr` (+ uniform 1/h optimisation)
  - [x] `lg_stencil` with endpoint reduction (numpy edge_order=1 default)
  - [x] option/spacing parse `lg_parse` (Method, DifferenceOrder, WindowLength, Axis, spacing)
  - [x] buffer kernel `list_gradient_ndarray` (float64/float32, per-axis strided axpy)
  - [x] List/symbolic path `lg_grad_axis_list` / `lg_grad_along_top`
  - [x] `builtin_list_gradient` entry + dispatch
  - [x] multi-axis assembly + Axis selection + numpy return-shape rules

## Registration
- [x] `src/list/list.h` — include list_gradient.h
- [x] `src/list/list_init.c` — add_builtin + ATTR_PROTECTED
- [x] `src/sym_names.{h,c}` — SYM_ListGradient declare/define/intern
- [x] `src/info.c` — docstring
- [x] `src/options_builtin.c` — Options[ListGradient] defaults

## Fast-path surfaces
- [x] `src/pack.c` — AWARE[] entry (NOT int64_ok)
- [x] `src/compile/compile_ndtables.c` — ND_FNS row + include

## Tests / docs
- [x] `tests/test_list_gradient.c` (24 cases) + `tests/CMakeLists.txt` (2 edits)
- [x] `docs/spec/builtins/arithmetic.md` (ListGradient section) + weekly changelog

## Verification
- [x] build clean; `list_gradient_tests` all pass; `differences_tests` regression pass
- [x] CompileDiagnostics array=True / scalar=False; compiled call + auto-compile Table
- [x] numpy parity spot-checked (1-D, non-uniform 17/6, 2-D, forward/backward, order4)
- [x] check-c99, check-packed-aware, check-nd-surfaces, check-array-exactness — clean
- [x] check-fastpath-sweep — ListGradient NOT flagged (float path fires). The
  audit's Error 1 is a PRE-EXISTING backlog of 25 unrelated heads (ChineseRemainder,
  Grad, Laplacian, EuclideanDistance, Xor, ...) not in OFF_BUFFER.
- [x] valgrind: byte-identical leak totals vs trivial script → zero ListGradient leaks
- note: check-compile-coverage and check-fastpath-sweep both have PRE-EXISTING
  backlogs of unrelated heads and are NOT CI gates (CI runs only check-c99,
  check-packed-aware, and the full glibc build — all green). ListGradient is not
  flagged by either.

## Review
- Added `ListGradient`, a numpy.gradient port, as one self-contained module.
  Reuses: `options_extract`/`OptEntry` (options), `nd_rotate_axes`/`ndred_accumulate`
  stride patterns (buffer walk), `expr_new_ndarray_like` (representation-preserving
  result), `ndstruct_delist_repack` (declined-NDArray → exact/symbolic bridge), and
  the `ND_FNS` table (Compile). Single Fornberg kernel underlies every
  scheme/order/window/grid; two instantiations (double + Expr) serve the machine
  and exact/symbolic paths.
