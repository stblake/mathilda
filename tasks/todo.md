# Task: Implement `AlgebraicNumberTrace`

Additive sibling of `AlgebraicNumberNorm`. Trace = sum of roots of minimal
polynomial = `-c_{n-1}/c_n`; relative trace over `Q(theta)` = `(n/d)·absolute`.

## Implementation
- [ ] Engine: `qqbar_abs_trace` + `flint_qqbar_algebraic_number_trace` + `#else` stub in `src/poly/flint_qqbar.c`
- [ ] Prototype + doc comment in `src/poly/flint_qqbar.h`
- [ ] New wrapper `src/poly/algebraicnumbertrace.c` + `.h`
- [ ] Register in `src/core.c` (fwd-decl + call)
- [ ] `SYM_AlgebraicNumberTrace` in `src/sym_names.h` + `.c`
- [ ] `Extension -> None` default in `src/options_builtin.c`
- [ ] Docstring in `src/info.c`
- [ ] Bump `src/version.h` 0.159 -> 0.160

## Tests
- [ ] New `tests/test_algebraicnumbertrace.c` (mirror norm test)
- [ ] CMake: COMMON_SRC + test target block in `tests/CMakeLists.txt`

## Docs
- [ ] `docs/spec/builtins/algebra.md` new section
- [ ] `docs/spec/changelog/2026-09-21.md` append entry

## Verification
- [ ] `make -j` + REPL smoke test of all spec examples
- [ ] Build + run unit tests (ctest)
- [ ] valgrind clean
- [ ] `make check-c99`, `make check-fastpath-sweep`
- [ ] Rebuild code-review-graph
- [ ] Commit + tag v0.160

## Review

Implemented `AlgebraicNumberTrace` as a precise mirror of `AlgebraicNumberNorm`.

- **Engine** (`src/poly/flint_qqbar.c`): `qqbar_abs_trace` reads coeff `x^{n-1}`
  and the leading coeff, giving `-c_{n-1}/c_n` (sum of roots); the relative case
  scales the absolute trace by the tower index `n/d` with `fmpq_mul_si` (norm
  raises to that power). Reuses `to_qqbar`, `algint_generator`,
  `qqbar_express_in_field`, `expr_from_fmpq`. `#else` stub added for `USE_FLINT=0`.
- **Wiring**: wrapper `algebraicnumbertrace.c/.h`, `core.c` init, `sym_names.{c,h}`,
  `Extension -> None` default in `options_builtin.c`, docstring in `info.c`,
  version bumped 0.159 → 0.160.

### Verification results
- Build: clean under `gcc-16 -std=c99 -Wall -Wextra` with FLINT/MPFR/LAPACK/etc.
- REPL smoke test: **all 12 spec examples produce the exact expected output**
  (10, -11, 2, -2/3, 0, -1, 2, 1, {10,0}, 10, -3, True), plus both error messages.
- Unit tests: `test_algebraicnumbertrace` all passed; 6 sibling algebraic-number
  test suites still pass (no regression).
- `make check-c99`: rc=0 (clean).
- valgrind: leak totals **byte-identical** to the proven-clean `AlgebraicNumberNorm`
  on the same inputs (13,440/420 def, 6,312/60 indir) → zero new leaks; remainder
  is the documented qqbar/FLINT/startup baseline noise.
- No leak stack trace mentions any AlgebraicNumberTrace symbol.
- code-review-graph rebuilt.
- `make check-fastpath-sweep`: the gate DID flag Trace (and, it turned out, Norm —
  a silent omission in the prior commit that left the gate red). Both are algebraic-
  number heads whose per-element cost is a qqbar min-poly read, not a numeric buffer
  op, so both were added to `SKIP_EXPLOSIVE` (the documented sibling convention,
  matching Denominator/Polynomial/IntegralBasis). Scoped `--only` check confirms both
  are now skipped (exit 0). Full gate re-run: neither Trace nor Norm appears in the
  NEW list — my change adds zero red heads and removes one (the overlooked Norm). The
  gate still exits 1 from a PRE-EXISTING 28-head backlog (LinearModelFit, NMaximize,
  Predict, ToNumberField, AlgebraicIntegerQ, Xor, ... from a stale OFF_BUFFER last
  re-recorded 2026-08-03) — unrelated to this change and out of scope.
