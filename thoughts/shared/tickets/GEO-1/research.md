---
ticket: GEO-1
created: 2026-08-27T17:41:03-0400
researcher: msollami
source_sha: 15b088dab072a89d29e2e8979f4e266de7cc12dc
branch: main
repository: mathilda
topic: "Add core computational-geometry builtins (Area, Perimeter, RegionCentroid, RegionMember on Polygon; ConvexHullRegion) in the spirit of WL 12 core geometry"
tags: [research, codebase, geometry, polygon, region, convex-hull, builtins, tests]
subsystems: [geometry, graphics, tests]
type: research
lifecycle: active
status: complete
last_updated: 2026-08-27
last_updated_by: msollami
---

# Research: Core computational-geometry builtins for Mathilda

**Date**: 2026-08-27T17:41:03-0400
**Researcher**: msollami
**Git Commit**: 15b088dab072a89d29e2e8979f4e266de7cc12dc
**Branch**: main
**Repository**: mathilda

## TL;DR
*(Human, ≤80 words)*
Asked: what exists for computational geometry, and what does adding WL-12-style region builtins require? Found: nothing user-facing exists — `Polygon` and friends are inert rendering tags; the only shoelace code is a private winding test in the renderer. All infrastructure needed (builtin registration, exact rationals via GMP, symbolic `Sqrt` sums, point-list loaders, test harness) exists and is well-trodden. Points toward: a new `src/geometry.c` module with five builtins. Unresolved: nothing blocking; scope decisions recorded as assumed (beta test, no human present).

## Summary
*(Human, ≤200 words)*
Mathilda has zero computational-geometry evaluation today. `Polygon`, `Point`, `Line`, `Circle`, `Disk`, `Rectangle` are interned symbols consumed only by `draw_primitive()` in the renderer; `Area`/`Perimeter`/`RegionCentroid`/`RegionMember`/`ConvexHullRegion` appear nowhere in `src/sym_names.h` (complete inventory of names the C code tests for). Baseline binary confirms all five heads return unevaluated.

Everything needed to add them exists: `symtab_add_builtin` with a NULL-decline ownership contract; module `_init()` wiring in `core_init()`; the makefile wildcard picks up any new `src/*.c` automatically; `make_rational`/`make_rational_mpz` build canonical exact results; the evaluator simplifies symbolic `Plus`/`Sqrt` (verified: `1+1+Sqrt[2]` → `2 + Sqrt[2]`, `Sqrt[8]` → `2 Sqrt[2]`), enabling WL-faithful exact `Perimeter`; `na_load_matrix` handles both nested `List` and visible `NDArray` point matrices. Packed lists are materialized before any non-AWARE head runs, so a new geometry head sees boxed exprs unless we opt in.

Reference semantics were pinned against a live Wolfram kernel (12 probe cases, exact and machine). The recommended slice: `src/geometry.c` with `Area`, `Perimeter`, `RegionCentroid`, `RegionMember` (2D `Polygon`) and `ConvexHullRegion` (returns `Polygon`), exact via GMP `mpq_t` with a machine-double path.

## Open Questions
*(Human, no cap — a question list, not prose)*

### Unresolved
_None._

### Resolved
- [x] Which vertical slice? — Five builtins: `Area`, `Perimeter`, `RegionCentroid`, `RegionMember` on 2D `Polygon[{{x,y},...}]`, plus `ConvexHullRegion[{pts}]` → `Polygon[hull]`. `ConvexHullMesh`/`MeshRegion` machinery is out of scope — Mathilda has no mesh-object type and WL 12.1's `ConvexHullRegion` returning a `Polygon` for 2D input is verified real behavior. (assumed — beta test)
- [x] Exact or machine arithmetic? — Both, WL-faithfully: if every coordinate is Integer/BigInt/Rational, compute exactly with GMP `mpq_t` (canonicalized via `make_rational_mpz`); any Real coordinate switches the whole computation to machine doubles. Verified against WL: `Area[Polygon[{{0,0},{1,0},{1/2,1/2}}]]` → `1/4` exact; `Area[Polygon[{{0,0},{1.5,0},{1.5,1},{0,1}}]]` → `1.5`. (assumed — beta test)
- [x] Self-intersecting polygons? — Out of scope; shoelace semantics apply (documented limitation). WL computes the enclosed-region area, which differs; simple (non-self-intersecting) polygons only for GEO-1. (assumed — beta test)
- [x] 3D polygons / higher-dim regions? — Out of scope for GEO-1; 2D only, decline (return NULL) on anything else so expressions stay unevaluated rather than wrong. (assumed — beta test)
- [x] Packed/NDArray/Compile surfaces (CLAUDE.md hard rule)? — Geometry heads are not element-wise numeric maps: they consume a `Polygon`-wrapped structural argument and return a scalar/structure. Packed lists are auto-materialized before non-AWARE heads run (`src/pack.c:472-478` no-nesting invariant), so correctness is safe by default. A visible `NDArray` inside `Polygon[...]`/`ConvexHullRegion[...]` will be read via the `is_ndarray` + `na_load_matrix` path so it is never left unevaluated. If `make check-packed-aware`/`check-compile-coverage` audits flag the new heads, add EXEMPT entries with one-line reasons — the documented, deliberate mechanism per CLAUDE.md. (assumed — beta test)
- [x] Does any existing code compute polygon area? — Yes, `polygon_signed_area()` (`src/graphics/render.c:90-97`, decl `src/graphics/render.h:85-93`), but it is double-only, renderer-private, and used solely for winding detection. New module computes its own exact/machine versions; no reuse (would invert the dependency: core module depending on optional `USE_GRAPHICS` code).

## Requires Approval
*(Human, ≤100 words)*
No human present in this beta run; all scope calls above are recorded "(assumed — beta test)" and a reviewer should confirm: (1) the five-builtin slice, (2) simple-polygon-only limitation, (3) `ConvexHullRegion` returning a plain `Polygon[...]` (WL attaches a cell spec `{1,2,...}` second argument; we return the bare vertex form), (4) `Area[Polygon[<2 pts>]]` → `Undefined` matching WL.

---

## Research Question
"Add core computational-geometry builtins to mathilda, in the spirit of Wolfram Language 12's core geometry (ConvexHullMesh-style hull computation, RegionMember-style point-in-region tests, area/perimeter/centroid of polygons — pick the smallest coherent vertical slice your research supports). Determine what geometry support already exists before assuming greenfield."

## Detailed Findings

### Baseline (reproduced, verbatim)
Build: `make -j8` with `SDKROOT=$(xcrun --show-sdk-path)`, CC auto-detected `gcc-16` → exit 0, `./Mathilda` produced (commit 15b088da).
`./Mathilda -file baseline_geom.m` output — all heads inert:
```
Area[Polygon[{{0, 0}, {2, 0}, {2, 1}, {0, 1}}]]
Perimeter[Polygon[{{0, 0}, {1, 0}, {0, 1}}]]
RegionCentroid[Polygon[{{0, 0}, {1, 0}, {0, 1}}]]
RegionMember[Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}], {1, 1}]
ConvexHullRegion[{{0, 0}, {1, 0}, {0, 1}, {2, 2}}]
```
Capability probes: `1/4` prints `1/4` (`FullForm` `Rational[1, 4]`); `Sqrt[2] + 2` → `2 + Sqrt[2]`; `Sqrt[8]` → `2 Sqrt[2]`; `Undefined` is a printable symbol; `N[Sqrt[2]+2]` → `3.41421`.

### Existing geometry support (none user-facing)
- `SYM_Polygon` interned at `src/sym_names.c:723,1605` (decl `src/sym_names.h:735`); consumed in exactly one place: `draw_primitive()` `src/graphics/render.c:1078-1109` (raylib `DrawTriangleFan`).
- `polygon_signed_area()` — shoelace, double-only — `src/graphics/render.c:90-97`, used only for winding direction at `render.c:1088-1103`.
- `Point`/`Line`/`Rectangle`/`Circle`/`Disk` interned `src/sym_names.c:716-722`/`1600-1604`; all pure drawing tags in `draw_primitive()` (`render.c:1019-1077`).
- No `Area`, `Perimeter`, `Centroid`, `RegionCentroid`, `RegionMember`, `ConvexHull*`, `Region` in `src/sym_names.h` (full 970-line inventory read; the file's own header states every strcmp'd name lives there). `SYM_RegionFunction` (`sym_names.h:812`) is a plot option, unrelated.
- `src/imagegeom.c` is raster resampling (`ImageResize`), no vector geometry.
- No geometry category in `Mathilda_spec.md` index (lines 27-62) nor `docs/spec/builtins/`.

### Builtin registration pattern (the recipe to follow)
- `typedef Expr* (*BuiltinFunc)(Expr* res);` `src/symtab.h:40`; `void symtab_add_builtin(const char*, BuiltinFunc);` `src/symtab.h:175`.
- Ownership contract (`SPEC.md` §4, `docs/extending.md:22-29`): return `NULL` to decline WITHOUT freeing `res` (expression stays unevaluated); return a fresh `Expr*` on success (evaluator frees `res`). Never `expr_free(res)` inside the builtin.
- Minimal module exemplar: `contfrac_init()` `src/contfrac.c:870-877` — `symtab_add_builtin(...)` then `symtab_get_def("Name")->attributes |= ATTR_...;`. Docstring-carrying exemplar: `ml_init` `src/ml/pca.c:429-463` with `symtab_set_docstring` (`src/symtab.h:178`).
- Wiring: `#include "geometry.h"` in `src/core.c` (~line 127 block) + `geometry_init();` call inside `core_init()` (`src/core.c:225`, sibling calls at 244-246). `sym_names_init()` runs first (`core.c:234`).
- New SYM_ names: decl in `src/sym_names.h`, `NULL` definition + `intern_symbol` in `sym_names.c` (`sym_names.c:12-13`, `:902+`); per repo CLAUDE.md, internal symbols must be defined there.
- Makefile: `SRC = $(wildcard src/*.c) ...` (`makefile:361-362`) — a new `src/geometry.c` is picked up automatically.
- Head checks: `head_is(e, SYM_List)` `src/common.h:38-52`; pointer-equality on interned names.
- Docs contract (repo CLAUDE.md): docstring via `symtab_set_docstring`, spec entry under `docs/spec/builtins/`, weekly changelog `docs/spec/changelog/<Monday>.md`, attributes assigned.

### Numeric extraction and result construction
- `na_read_scalar`/`na_load_vector`/`na_load_matrix` (`src/linalg/numarray.h:42-72`, impl `numarray.c:26-282`): accept visible `NDArray` OR rectangular nested `List`, handle Integer/BigInt/Real/MPFR/`Rational[p,q]`/`Complex`, reject ragged input. Right tool for the machine path and for NDArray acceptance.
- Exact path: `make_rational(int64_t,int64_t)` `src/arithmetic.c:34-53`; `make_rational_mpz(mpz_t,mpz_t)` `arithmetic.c:61-84` (canonicalizes, demotes to Integer when denominator 1). `expr_to_mpz`, `expr_is_integer_like`, `expr_is_numeric_like` (`src/expr.h:336-340`). GMP `mpq_t` available (gmp linked unconditionally).
- Constructors: `expr_new_integer`/`expr_new_real`/`expr_new_bigint_from_mpz` (`src/expr.h:226-268`); lists via `expr_new_function(expr_new_symbol(SYM_List), args, n)` (exemplar `src/contfrac.c:106-113`); symbolic results built as expression trees then `eval_and_free` (exemplar `fcf_qirr_to_expr` `src/contfrac.c:722-756` builds `Plus`/`Times`/`Sqrt`/`Power`). `True`/`False` = interned `SYM_True`/`SYM_False`; predicates must never return NULL (`docs/extending.md:85-86`).
- Point-list-consuming builtin exemplars: `PrincipalComponents` (`src/ml/pca.c:327-364`, via `ml_read_data` → `na_load_matrix`), `Fit` (`src/fit.c`, own loader + local `fit_to_double` `src/fit.c:138-159`).

### Packed-array / NDArray / Compile surfaces
- `AWARE` list `src/pack.c:489-782`, `INT64_OK` `pack.c:873+`. No graphics/geometry head on either.
- Transparency gate: packed (invisible) List materialized into boxed Exprs before any non-AWARE head evaluates (`pack.h:38-39`); no-nesting invariant `pack.c:472-478` means a packed buffer can never sit nested inside `Polygon[...]`'s boxed tree.
- Visible `NDArray[...]` is NOT intercepted: `ConvexHullRegion[NDArray[...]]` or `Polygon[NDArray[...]]` would reach the builtin as a real `EXPR_NDARRAY` node → must be read explicitly (`is_ndarray` + `na_load_matrix`), as `Fit`/`ArrayPlot` do (`pack.c:701-719,774-781`).
- Audits: `make check-packed-aware` → `tools/check_packed_aware.py` (EXEMPT dict at :47 with per-entry reasons); `check-compile-coverage` → `tools/compile_coverage.py` (EXEMPT :186, BASELINE ratchet :263, mirrored in `COMPILE_MISSING.md`); `check-array-exactness` → `tools/check_array_exactness.py` (EXEMPT with Mathematica-output justification; needs built binary, ~5 min); `check-nd-surfaces`, `check-fastpath-sweep` (ratchet `OFF_BUFFER` at `tools/nd_fastpath_sweep.py:608`).

### Test harness
- `tests/CMakeLists.txt`: kernel sources compiled once as OBJECT library `mathilda_common` from `COMMON_SRC` (entries lines 281-887; `add_library` at :899). **A new `src/geometry.c` must be added to `COMMON_SRC`** or every test binary fails to link (the known trap). New test binary block exemplar (`:942-945`):
  ```cmake
  add_executable(trigexp_zero_tests test_trigexp_zero.c $<TARGET_OBJECTS:mathilda_common>)
  target_link_libraries(trigexp_zero_tests m)
  target_include_directories(trigexp_zero_tests PRIVATE ../src)
  add_test(NAME trigexp_zero_tests COMMAND trigexp_zero_tests)
  ```
  `add_test` is mandatory for ctest visibility (`enable_testing()` at :31; several legacy targets lack it and are invisible).
- Harness: `assert_eval_eq(const char* input, const char* expected, int is_fullform)` `tests/test_utils.h:18-30` (string-compare on printed result); `assert_eval_startswith` :35-50; `ASSERT`/`ASSERT_STR_EQ`/`ASSERT_MSG` (NDEBUG-safe, exit(1)); `TEST(name)` macro :52; 60s alarm watchdog :101-106. Setup: `symtab_init(); core_init();` in `main()`; no teardown. Exemplar file: `tests/test_trigexp_zero.c` (151 lines); also memory-loop tests for valgrind coverage (`tests/test_clip.c` pattern, per 2026-08-18 UnitBox research).
- Run: `cd tests && mkdir -p build && cd build && cmake .. && make -j8 && ctest` (or run `./geometry_tests` directly). GCC-only enforced in `tests/CMakeLists.txt:1-25` too.
- Non-interactive scripts: `-file` handled `src/repl.c:985-1009` → `run_script_file` :910-963; output only via `Print[]`; exit 0/1. Piped stdin is NDJSON mode, not REPL.

### Verification gates (for the ladder)
- `make check-c99` → `tools/check_c99_portability.py` (no exemptions; new POSIX symbol = one line in its FUNCTIONS table). CI runs it before build.
- CI (`.github/workflows/build.yml`): check-c99 → check-packed-aware → `make -j` → NDJSON stdin smoke test; second job builds with all optional libs off (`USE_MPFR=0 USE_FLINT=0 USE_ECM=0 USE_LAPACK=0 USE_GRAPHICS=0`) — geometry.c must not depend on optional libs (GMP/readline are unconditional; raylib/graphics is NOT).
- No `make test` target; ctest is the unit rung. No general linter; `-Wall -Wextra -Werror=…` promotions are the static rung (makefile CFLAGS).

### Reference semantics (live Wolfram kernel, verified 2026-08-27)
```
Area[Polygon[{{0,0},{1,0},{1/2,1/2}}]]                      -> 1/4
Perimeter[Polygon[{{0,0},{1,0},{0,1}}]]                     -> 2 + Sqrt[2]
RegionCentroid[Polygon[{{0,0},{1,0},{0,1}}]]                -> {1/3, 1/3}
RegionMember[Polygon[{{0,0},{2,0},{2,2},{0,2}}], {2,1}]     -> True   (boundary inclusive)
RegionMember[...same square...], {1/2,1/2}]                 -> True
ConvexHullRegion[{{0,0},{1,0},{0,1},{2,2},{1/2,1/2}}]       -> Polygon[{{0,0},{1,0},{2,2},{0,1}}, {1,2,3,4}]
Area[Polygon[{{0,0},{1.5,0},{1.5,1},{0,1}}]]                -> 1.5
Area[Polygon[{{0,0},{4,0},{4,4},{2,1},{0,4}}]]              -> 10        (concave, exact)
RegionCentroid[Polygon[{{0,0},{4,0},{4,4},{2,1},{0,4}}]]    -> {2, 7/5}
RegionMember[Polygon[{{0,0},{4,0},{4,4},{2,1},{0,4}}], {2,3}] -> False  (in notch)
ConvexHullRegion[{{0,0},{2,0},{1,0},{2,2},{0,2},{1,1}}]     -> Polygon[{{0,0},{2,0},{2,2},{0,2}}]  (dedups collinear/duplicate/interior)
Perimeter[Polygon[{{0,0},{3.,0},{3.,4.}}]]                  -> 12.
Area[Polygon[{{0,0},{1,0}}]]                                -> Undefined
```
Note: WL's `ConvexHullRegion` output carries a second cell-spec argument `{1,2,...,n}`; GEO-1 returns the bare `Polygon[{verts}]` form (recorded under Requires Approval).

## Code References
- `src/symtab.h:40,175,178` — BuiltinFunc type, symtab_add_builtin, symtab_set_docstring
- `src/core.c:225,234,244-246` — core_init, sym_names_init first, sibling module init calls
- `src/contfrac.c:870-877` — minimal module init exemplar; `:550-605` decline-with-NULL exemplar; `:722-756` symbolic-result construction
- `src/sym_names.c:716-724,1598-1606` — graphics-primitive interning site (where new SYM_ entries go alongside)
- `src/arithmetic.c:34-53,61-84` — make_rational / make_rational_mpz
- `src/linalg/numarray.h:42-72`, `numarray.c:26-282` — NDArray/List matrix loader
- `src/common.h:38-52` — head_is
- `src/graphics/render.c:90-97,1078-1109` — existing shoelace (private) and Polygon rendering
- `src/pack.c:472-478,489-782,701-723` — no-nesting invariant, AWARE list, NDArray-visible-argument precedent
- `tests/CMakeLists.txt:281-887,899,942-945` — COMMON_SRC, object library, test block exemplar
- `tests/test_utils.h:18-57,101-106` — assert_eval_eq and NDEBUG-safe asserts, watchdog
- `src/repl.c:985-1009,910-963` — -file script execution
- `makefile:361-366,469-612` — SRC wildcard, check-* targets
- `.github/workflows/build.yml:26-158` — CI order, no-optional-libs job

## Architecture Insights
- Graphics primitives are deliberately inert; giving `Area`/`RegionMember` their own builtins (rather than making `Polygon` evaluate) preserves that design — `Polygon` stays a structural tag, geometry heads pattern-match on it, exactly as `draw_primitive` does.
- Exact-vs-machine duality is a repo-wide idiom (three separate leaf→double families exist deliberately unshared); a geometry module should follow the local-helper convention rather than adding cross-module deps on `numarray` for the exact path (but SHOULD use `na_load_matrix` for the NDArray machine path, as `Fit` does).
- Core modules must not depend on `USE_GRAPHICS`-gated code (CI builds with it off) — hence no reuse of `render.c`'s shoelace.
- Predicates (`RegionMember`) must always return `True`/`False` for well-formed input, never decline into an unevaluated form mid-computation, matching the `*Q` convention — but malformed/symbolic input should decline (return NULL) like every other builtin.

## Historical Context (from thoughts/)
- `thoughts/shared/research/2026-08-18-unitbox-builtin.md` — the closest prior "add a builtin" research: confirms the piecewise.c registration/attribute pattern, `tests/test_clip.c` as test exemplar (TEST macro + memory-loop for valgrind), and the CMake add_test gotcha (clip_tests built but not ctest-registered).
- `thoughts/shared/tickets/DEMO-2/` — prior ticket-folder-convention artifact (shape reference for this ticket).
- `thoughts/shared/research/2026-08-21-DEMO-1-nminimize-nmaximize-benchmarking.md` + `thoughts/shared/plans/2026-08-21-DEMO-1-nminimize-nmaximize-benchmarks.md` — prior flat-file-convention artifacts (not geometry-relevant).

## Related Research
- `thoughts/shared/research/2026-08-18-unitbox-builtin.md`
