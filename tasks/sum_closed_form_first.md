# Sum: prefer closed form over brute-force expansion on finite numeric ranges

## Problem

`Sum[3/2 i^2 - i/2, {i, -10000, 10000}]` takes ~0.16 s vs Mathematica's
~0.036 s. Root cause: `sum_one_spec` (src/sum/sum.c) unconditionally
brute-force expands any finite numeric range term-by-term (20001 evals +
a 20001-arg Orderless-sorted `Plus` fold), only consulting the existing
closed-form cascade when the span exceeds the 100M runaway guard.

The closed form already exists and is correct:
`Sum`Polynomial[3/2 i^2 - i/2, i, -10000, 10000]` -> 0.00035 s (~470x faster).
Secondary bug: `Method -> "Polynomial"` is ignored for numeric ranges
because expansion runs before method dispatch.

## Fix

In `sum_one_spec`, for a RANGE spec with a **unit positive integer step**
and **non-empty ascending bounds**, try the closed-form Method cascade
*before* falling back to `expand_range`. Telescoping
`F(imax+1) - F(imin)` is exact for any imin/imax differing by a
non-negative integer, so this is correctness-preserving; non-summable
bodies (cascade returns NULL) fall through to expansion as today.

Guards (all required for closed-form-first):
- `!is_real` (integer bounds/step) — avoids float edge cases.
- `di_val == 1` (unit step) — the unit-step antidifference is invalid for
  di != 1; those keep expanding.
- `min_val <= max_val` (non-empty) — empty ranges must yield 0 via
  expansion; the telescoping form would give a wrong nonzero value.

Must shadow the iterator var around the cascade call (Sum is HoldAll;
an outer binding of `i` would otherwise leak into the held body / stage
args). Expansion already shadows internally.

## Tasks

- [x] Diagnose (measured).
- [x] Edit `sum_one_spec` in `src/sum/sum.c`: closed-form-first block.
- [x] Build (`make -j8`), confirm clean.
- [x] Verify: original example fast + same answer; empty range, step != 1,
      non-polynomial finite range, symbolic bound, outer binding all correct.
- [x] Add 8 regression checks to `tests/test_sum.c`.
- [x] Run `sum_tests` binary only (scoped) — 42/42 pass.
- [x] Update `docs/spec/builtins/calculus.md` + 2026-05-25 changelog.

## Review

Single localized change in `sum_one_spec`: when a RANGE spec resolves to a
non-empty, integer-bounded, unit-step numeric span, try `dispatch_def`
(the existing Method cascade) first, shadowing the iterator var, and only
fall back to `expand_range` if the cascade returns NULL. No new files,
no API change, additive guards.

Measured: `Sum[3/2 i^2 - i/2, {i, -10000, 10000}]` 0.165 s -> 0.0003 s
(~500x), identical value 1000150005000. Gates verified: empty range ->0,
step 2 ->expansion, If[]/1/i/Prime[] bodies ->expansion, `i=7` outer
binding does not leak (->30). `sum_tests` 42/42.

Secondary fix: `Method -> "Polynomial"` now takes effect on finite
unit-step numeric ranges (was silently ignored before).

Follow-up (now fixed, 2026-05-30): the symbolic-bound dispatch path
dropped `s.di`, so a symbolic-bound sum with step != 1 (e.g.
`Sum[i, {i, 1, n, 2}]`) wrongly collapsed to the unit-step closed form
`1/2 n (1 + n)`. The closed-form stages take no step argument and assume a
unit step, so the final `dispatch_def` is now guarded by
`is_unit_step(s.di)`: a non-unit step with no step-aware closed form
returns NULL and the `Sum[...]` stays held. Two regression checks added
(`Sum[i, {i, 1, n, 2}]`, `Sum[i^2, {i, 1, n, 3}]`); `sum_tests` 44/44.
