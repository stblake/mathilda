# LatticeReduce (LLL lattice basis reduction)

## Goal
`LatticeReduce[{v1, v2, ...}]` returns an LLL-reduced basis for the lattice
spanned by the row vectors. Entries may be integers, Gaussian integers,
rationals, or Gaussian rationals. Exact arithmetic (GMP `mpq_t`), so it works
for machine-size and arbitrary-precision (bignum) inputs alike.

## Algorithm
Classic LLL (Lenstra–Lenstra–Lovász), delta = 3/4, with **exact rational
Gram–Schmidt** generalized to the complex (Hermitian) case so one code path
covers all four number kinds:
- Gaussian rational entry = `{mpq_t re; mpq_t im}` (`GRat`).
- Inner product `<x,y> = sum x_k conj(y_k)`.
- Size reduction rounds `mu[k][l]` to the nearest Gaussian integer.
- Incremental GSO: Gram–Schmidt once; size-reduction updates `mu` in place;
  swaps update `mu`/`Bnorm` via the conjugate-aware Cohen swap formulas
  (no full recompute). Efficient.
- Lattice preserved exactly (all ops are Z/Z[i] row combos + swaps), so
  `Abs[Det]` and right null-space relations are invariant.
- Independent rows required (every documented use case is full rank);
  dependence detected at GSO time -> diagnostic + unevaluated.

## Steps
- [x] Explore linalg / arithmetic / test conventions.
- [ ] `src/linalg/latticereduce.c` — GRat helpers, conversions, LLL, builtin.
- [ ] `src/linalg/linalg.h` — declare `builtin_latticereduce`.
- [ ] `src/linalg/linalg.c` — register builtin + ATTR_PROTECTED.
- [ ] `src/info.c` — terse docstring.
- [ ] `src/sym_names.{h,c}` — intern `LatticeReduce`.
- [ ] `tests/test_latticereduce.c` + `tests/CMakeLists.txt`.
- [ ] Build, run scoped test, valgrind for leaks.
- [ ] Docs: `docs/spec/builtins/...` + weekly changelog.

## Validation / errors (per spec)
- `argx`   — arg_count != 1.
- `matrix` — not a non-empty rectangular matrix.
- `latm`   — an entry is not rational / Gaussian-rational.

## Review (2026-06-04 — complete)

All steps done. `LatticeReduce` reproduces the documented Wolfram output
**exactly** on every reference example — including the 10⁸ and 10²⁰ bignum
integer-relation lattices (`{{-3,0,0,1,0},...}` and first row `{0,1,-4,1,0}`),
the 2-D `{{12,2},{13,4}} -> {{1,2},{9,-4}}`, the 3×4 relation matrix, and the
3×4 `1345/35/154` lattice `{{0,9,-2,7},{1,1,-9,-6},{1,-3,-8,8}}`.

- Exact GMP-rational LLL (δ = 3/4); no floating point → correct for machine
  and bignum inputs alike. Hermitian Gram–Schmidt → integers, rationals,
  Gaussian integers, Gaussian rationals share one path.
- Incremental GSO (in-place size reduction + conjugate-aware Cohen swap) —
  efficient, no per-swap recomputation.
- 19 unit tests pass (`tests/test_latticereduce.c`); `linalg_tests` regression
  clean.
- Valgrind: zero leak frames reference `latticereduce.c`; the residual
  "definitely lost" total matches the `matrank_tests` baseline (pre-existing
  harness/global allocations), so no leak is attributable to the new code.
- Files: `src/linalg/latticereduce.c` (new), `linalg.{c,h}`, `info.c`,
  `sym_names.{c,h}`, `tests/{test_latticereduce.c,CMakeLists.txt}`, docs spec
  + changelog.

### Known scope
Linearly independent rows required (every documented use case is full rank);
a dependent generating set is reported via `LatticeReduce::dep` rather than
returning zero-padded rows (MLLL) — a possible future extension.
