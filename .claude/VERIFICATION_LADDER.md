# Verification ladder configuration

Read by `skills/verification-ladder/scripts/ladder.py`. Pinned by `/ais:onboard` on
2026-08-31 from `detect_ladder.py`'s DETECT pass.

**The pinning rule used here: a rung takes a command only if DETECT could cite a real CI
line for it.** This repo's `detect_ladder.py` run produced candidates for four toolchains
— `c-cmake`, `make`, `node-ts`, `python` — but only two carried a `ci:` citation. Every
other candidate was a `hint:` derived from a manifest in a *subdirectory*
(`benchmarks/requirements.txt`, `frontend/package.json`, `tests/CMakeLists.txt`), so it
would have measured a different part of this repo than the C tree under `src/`. A rung
reporting green because it linted `benchmarks/` is worse than a rung honestly reporting a
gap, so the hints are rejected and recorded below rather than silently dropped.

`not-configured` here is a confirmed "nothing exists for this rung yet", not a
placeholder. `ladder.py` reports it as a coverage gap, never a pass — and because this
file exists, it no longer falls back to probing for a tool on `PATH`.

## Static analysis

```
static = not-configured
static-reviewed-against = clang-tidy -p build $(git ls-files '*.c' '*.cc' '*.cpp' '*.h' '*.hpp') | cppcheck --error-exitcode=1 . | eslint . | ruff check . | flake8 .
```

Rejected, all five: `clang-tidy` and `cppcheck` were `hint:c-cmake` from
`tests/CMakeLists.txt` (subdirectory-only, not this project's primary build system — the
root build is the `makefile`); `eslint` was `hint:node-ts` from `frontend/package.json`;
`ruff` and `flake8` were `hint:python` from `benchmarks/requirements.txt`. No CI citation
for any of them. Nothing proposed in their place — see the pinning rule above.

## Type checking

```
typecheck = not-configured
typecheck-reviewed-against = cmake --build build | tsc --noEmit | mypy . | pyright
```

Rejected, all four, same reason — every one a `hint:` from a subdirectory manifest, none
CI-cited. Note that for C this rung is normally "does it build", and this repo *does* have
a real build (`make -j"$(nproc)"`, `.github/workflows/build.yml`) — but `detect_ladder.py`
did not surface it as a typecheck candidate, so pinning it here would be an invention
rather than a confirmation. Left as an honest gap; a later pass can add it deliberately.

## Unit tests

```
unit = make check-c99 && make check-packed-aware
unit-reviewed-against = ctest --test-dir build --output-on-failure | make check-c99 | make check-packed-aware | python3 -m pytest -q -m "not integration" | pytest -q
```

The only two accepted commands in this file, and the only two DETECT could cite:

- `make check-c99` — `ci:.github/workflows/build.yml:73` ("Portability gate")
- `make check-packed-aware` — `ci:.github/workflows/build.yml:78` ("Packed-array opt-in
  audit")

Both run *before* the `Build` step in that workflow, and the workflow's own comment states
`check-packed-aware` is "source-level, so it needs no build". So neither needs a built
`./Mathilda` and this rung is runnable from a clean tree.

Rejected from this rung: `ctest --test-dir build` (`hint:c-cmake`, subdirectory-only, and
would need a configured CMake build dir) and both `pytest` invocations
(`hint:python`, from `benchmarks/`).

One honest caveat about placement: both accepted commands are *source-level audits*, which
is static-analysis-shaped work sitting on the `unit` rung. They are pinned here because
that is where DETECT proposed them and where they were confirmed. The C unit suite
(`tests/`, built with CMake) is genuinely not covered by this rung.

## Integration tests

```
integration = not-configured
integration-reviewed-against = python3 -m pytest -q -m integration
```

Rejected: `hint:python` from `benchmarks/requirements.txt`. This repo has no integration
tier for the C tree, so `not-configured` is the correct answer rather than a gap to close.
