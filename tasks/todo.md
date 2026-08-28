# ArrayPad + ArrayReshape implementation

Plan: `/Users/user/.claude/plans/let-s-implement-arrayreshape-in-luminous-pine.md`

## Tasks

- [x] 1. Shared 1D scheme helper `src/list/pad_schemes.{c,h}` (all 8 named schemes + constant/cyclic)
- [x] 2. `src/list/array_pad.{c,h}` — builtin_array_pad (amounts parse, options, recursion, mindimsize)
- [x] 3. `src/list/array_reshape.{c,h}` — builtin_array_reshape
- [x] 4. `src/ndstruct.{c,h}` — ndstruct_arraypad + ndstruct_arrayreshape (rank-1 fast paths)
- [x] 5. Registration: list.h, list_init.c, attr.c, info.c, sym_names.{h,c}
- [x] 6. ND registration: pack.c AWARE + INT64_OK
- [x] 7. Compile ratchet: tools/compile_coverage.py BASELINE
- [x] 8. Docs: lists-and-iteration.md + changelog 2026-08-24.md
- [x] 9. Tests: test_array_pad.c, test_array_reshape.c + CMake wiring
- [x] 10. Verify: build, check-c99, unit tests, ND audits, leaks

## Review

Implemented `ArrayReshape` and `ArrayPad` completely (all amount forms + all 8
named padding schemes), sharing a 1-D scheme engine (`pad_schemes.c`).

Key design points:
- Two builders in `array_pad.c`: a constant/cyclic recursive builder (fills
  full-width padding at every level, using `orig_dim + lo + hi` targets like
  `pr_build`), and a value-dependent outer-to-inner builder (extends each axis'
  fiber by the scheme, then recurses). `Extrapolated` uses the outward binomial
  recurrence and must `Expand` its combos to match Wolfram's expanded
  polynomials; the antisymmetric `*Differences` schemes fold `2*edge - x` at each
  boundary reflection.
- Packed fast paths (`ndstruct.c`) decline any non-exact shape/fill to the List
  path; default integer `0` into a real array is a legitimate mixed result
  (Wolfram does the same), recorded EXEMPT in `check_array_exactness.py`.

Verification: both new unit suites pass; `make check-c99`, `check-packed-aware`,
`check-array-exactness` green; sibling suites (`constant_array`, `array_flatten`,
`ndarray`) unregressed; 0 leaks under macOS `leaks`. Every reference output from
the user's spec reproduced exactly.

Pre-existing, unrelated: `make check-compile-coverage` is red on ~22 `Image*` /
`Interpolation` heads from earlier graphics commits that lack BASELINE entries;
`ArrayPad`/`ArrayReshape` are correctly baselined and not in that failure set.
