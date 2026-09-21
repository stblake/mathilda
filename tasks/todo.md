# Task: Implement `AlgebraicNumberPolynomial[a, x]`

Plan: `/Users/user/.claude/plans/continuing-on-from-our-crystalline-ladybug.md`

## Steps

- [x] Create `src/poly/algebraicnumberpolynomial.h`
- [x] Create `src/poly/algebraicnumberpolynomial.c` (structural, no FLINT)
- [x] Wire `algebraicnumberpolynomial_init()` into `src/core.c`
- [x] Add `SYM_AlgebraicNumberPolynomial` to `src/sym_names.h` + `src/sym_names.c`
- [x] Add docstring in `src/info.c` (terse, no examples)
- [x] Add attributes Listable | Protected (in _init)
- [x] Create `tests/test_algebraicnumberpolynomial.c`
- [x] Register test in `tests/CMakeLists.txt` (COMMON_SRC + add_test)
- [x] Create `tests/scripts/algebraicnumberpolynomial_leakcheck.sh`
- [x] Add `AlgebraicNumberPolynomial` to `SKIP_EXPLOSIVE` in `tools/nd_fastpath_sweep.py`
- [x] Docs: `docs/spec/builtins/algebra.md` entry
- [x] Docs: `docs/spec/changelog/2026-09-21.md` (created) + row in `Mathilda_spec.md` table
- [x] Bump `src/version.h` 0.155 -> 0.156
- [x] Build clean, run REPL examples, run unit tests, leak check, audits
- [x] Rebuild code-review graph

## Review

Implemented `AlgebraicNumberPolynomial[a, x]` — the inverse of `AlgebraicNumber`:
for `a = AlgebraicNumber[theta, {c0,...,cn}]` it returns `c0 + c1 x + ... + cn x^n`.

**Design:** purely structural (reads the coefficient vector already stored in the
object and builds a `Plus`/`Times`/`Power` tree, canonicalised by the evaluator).
No FLINT is used — it provides no benefit here, and the head therefore works in a
`USE_FLINT=0` build. Integer/rational inputs pass through unchanged; other inputs
emit `AlgebraicNumberPolynomial::naobj` and stay unevaluated. Attributes
`{Listable, Protected}`.

**Verification (all green):**
- Every example from the request reproduced in the REPL, incl. the round-trip
  `b = poly /. x -> Sqrt[2+Sqrt[3]]` → `1 + 2 Sqrt[2+Sqrt[3]] + 3(2+Sqrt[3]) +
  4(2+Sqrt[3])^(3/2)` and `RootReduce[b == a]` → `True`, and the addition-via-
  polynomials `2 + 4 x^2 + 5 x^3`.
- `tests/test_algebraicnumberpolynomial.c`: all passed (passthrough, build,
  Listable over both args, round-trip via RootReduce zero test, addition identity,
  naobj message via stderr capture, wrong-argc decline).
- Leak gate `algebraicnumberpolynomial_leakcheck.sh`: **0 leaks for 0 bytes**.
- `make check-packed-aware`: OK (head not flagged; no exemption needed — new file
  has no NDArray dispatch markers).
- `tools/nd_fastpath_sweep.py --only AlgebraicNumberPolynomial --gate-only`: 0
  shapes flagged (SKIP_EXPLOSIVE honored).
- `make check-c99`: clean.
- Regression: `algebraicnumber_tests`, `rootreduce_tests`, `minimalpolynomial_tests`
  all PASS.

**Not done (awaiting user):** git commit + tag `v0.156` (per CLAUDE.md, commit/tag
only on request).
