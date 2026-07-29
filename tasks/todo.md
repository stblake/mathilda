# Auto-compilation parity for the functional-programming heads

(Previous task — Limit oscillatory normal form — archived to
`plans/LIMIT_OSCILLATORY_TODO.md`.)

Commit `d42f543` brought the functional-programming heads into `Compile[]`'s
subset. This task closes the other half: the same heads called **directly**, with
no `Compile[]` wrapper, should take the automatic machine-number fast path
(`src/numloop.c`) the way `Nest`, `Map` and `Fold` already do.

## Audit — measured, `main` @ d1fb77b, n = 200000

`numloop` ON vs `MATHILDA_NO_NUMLOOP=1`, seconds:

| Head | Compile subset | Auto-compiled today | Time | Note |
|---|---|---|---|---|
| `Nest` | yes | **yes** (`numloop_nest`) | 0.0028 | 55x vs interp |
| `Map` | yes | **yes** (`numloop_map`) | 0.030 | 5.9x vs interp |
| `Fold` | yes | **yes** (`numloop_fold`) | 0.0076 | 12x vs interp |
| `FixedPoint` | yes | **yes** (`numloop_fixedpoint`) | 5e-6 | 5x vs interp |
| `NestWhile` | yes | **yes** (`numloop_nestwhile`) | 0.0027 | 46x vs interp |
| `Table` | yes | **yes** (`autocompile`, inexact iterator) | — | already wired |
| `NestList` | yes | **no** | 0.179 | 64x slower than `Nest` |
| `FoldList` | yes | **no** | 0.114 | 15x slower than `Fold` |
| `Scan` | yes | **no** | 0.138 | 4.6x slower than `Map` |
| `NestWhileList` | yes | **no** | 0.128 | — |
| `FixedPointList` | yes | **no** | (converges in 92) | no large case |
| `Accumulate` | yes | **no** | 0.067 | `evaluate(Plus[..])` per element |
| `Reverse`/`Sort`/`Take`/`Drop`/`Flatten`/`Transpose` | yes | n/a | 0.006–0.024 | no function body to compile; already native C, allocation-bound |

The five `*List` / `Scan` heads are the real gap: every hook site in `funcprog.c`
is guarded `if (!as_list)`, so the list-producing twin of each covered head falls
straight through to the interpreted loop.

## Plan

- [ ] 1. `numloop_nestlist(f, x0, n)` — n+1 Reals. Gate exactly as `numloop_nest`
      (machine-real seed, inexact-result guarantee, non-finite -> bail).
      Scalar seed only; an NDArray seed returns NULL (result would be a List of
      arrays, not the fused in-place map `numloop_nest_array` does).
- [ ] 2. `numloop_foldlist(f, x0, list)` — m+1 Reals, gated as `numloop_fold`.
- [ ] 3. `numloop_nestwhilelist(f, x0, test)` — {x0, ..., first value failing the
      test}, gated as `numloop_nestwhile`.
- [ ] 4. `numloop_fixedpointlist(f, x0)` — {x0, ..., fp, fp} (the final duplicate
      is in Mathilda's output: `FixedPointList[Cos,1.0]` has length 92 with
      `[[-1]] == [[-2]]`), gated as `numloop_fixedpoint`.
- [ ] 5. `numloop_scan(f, expr)` — run the compiled body per element for
      non-finite parity, return `Null`.
- [ ] 6. `numloop_accumulate(list)` — running sum in doubles, replacing the
      `evaluate(Plus[out[i-1], args[i]])` per element. Requires **every** element
      inexact, so the exact prefix of a mixed list cannot be flattened to a Real.
- [ ] 7. Wire each into `funcprog.c` / `accumulate.c`, dropping the `!as_list`
      guards.
- [ ] 8. Differential tests: fast path vs `MATHILDA_NO_NUMLOOP=1` must agree
      **bit-for-bit**, including the bail cases (non-finite, exact input, symbolic
      element, empty list).
- [ ] 9. `docs/spec/builtins/` + this week's changelog.

## Non-goals

- The structural heads (`Reverse`, `Sort`, `Take`, `Drop`, `Flatten`,
  `Transpose`) have no function argument. There is nothing to compile: they are
  already straight C over the argument list and their cost is the Expr
  allocation of the result. Their `Compile` support exists so they can appear
  *inside* a compiled body, where the values are unboxed machine arrays.
- MPFR / arbitrary-precision paths, as everywhere else in `numloop`.

## Findings noted along the way

- `numloop_nestwhile` and `numloop_fixedpoint` carry `CAP = 1000000`. On hitting
  it they discard **all** the compiled work and return NULL, so the interpreter
  redoes the whole loop: `NestWhile[#+1.&, 0., #<1000000.&]` costs 0.664s, versus
  0.0027s at 200k iterations. Pre-existing, correct-but-wasteful; recorded, not
  fixed here.

## Review

All nine items done. `src/numloop.c` gained `numloop_nestlist`,
`numloop_foldlist`, `numloop_nestwhilelist`, `numloop_fixedpointlist`,
`numloop_scan` and `numloop_accumulate`; each scalar twin became an
`*_impl(..., bool as_list)` with two thin wrappers, matching how `funcprog.c`
already spells `nest_impl(res, as_list)`. The two loops whose length is not known
until they run collect into a growable `DVec` and box to `Expr` only on success,
so a bail costs no allocation and has nothing to unwind.

Measured at n = 200000 (median of 3, `MATHILDA_NO_NUMLOOP=1` as the control):
NestList **6.6x**, Scan **11.1x**, NestWhileList **7.5x**, FoldList **5.0x**,
Accumulate **3.1x**.

Two bugs came out of it, neither one anticipated by the plan:

1. **The gate rounded exact values it was only passing through** — and the same
   hole was live on the shipped SCALAR paths, so `Map[# &, {1., 2, 3}]` answered
   `{1., 2., 3.}` against the interpreter's `{1., 2, 3}`, and `Fold[#2 &, 1.,
   {1, 2, 3}]` answered `3.` against `3`. Inexactness is a property of each
   result POSITION, not of the call: a computed value is Real if the body carries
   a Real literal, but a passed-through value keeps its input's type. The rule is
   now stated once in `numloop.h` and enforced per entry point.
2. **`VarCtx` was used uninitialised** in `compile_function` — three of seven
   fields set, and the optional `defined` pointer is dereferenced the moment it
   is non-NULL. Latent for as long as callers happened to leave zeros there;
   adding one more caller changed the stack shape and it faulted. Found by
   AddressSanitizer after `compile_tests` started segfaulting; `= {0}` fixes it.

Verified: a 3057-case differential corpus is byte-identical with the fast path on
and off; twelve test suites pass (numloop, compile, autocompile,
compile_coverage, compiledfunction, nestlist, foldlist, nestwhilelist,
fixedpointlist, list, funcprog_through, scan); `compile_tests` clean under ASAN;
`make check-c99` clean.

One test expectation of mine was wrong and worth recording: I asserted
`Module[{c = 0}, Scan[(c = c + 1) &, ...]; c] == 4`, which is Mathematica's
answer, not Mathilda's — a `Function` body does not close over a `Module`'s
renamed local here, so both paths answer 0. The test now uses a global counter
and measures the fast path rather than that (pre-existing) scoping behaviour.

Not done, deliberately: the `CAP = 1000000` waste in `numloop_nestwhile` /
`numloop_fixedpoint` noted above. It is correct, pre-existing, and orthogonal.
