---
source: src/core.c
---
**Algorithm.** `builtin_clip` accepts `Clip[x]` (clamp to `[-1, 1]`),
`Clip[x, {min, max}]`, and `Clip[x, {min, max}, {vmin, vmax}]`. A `List`
first argument is threaded manually — `Clip[{x1,...}, ...]` maps to
`{Clip[x1, ...], ...}` while the `{min,max}` / `{vmin,vmax}` configuration lists
are carried through unchanged (this is the Listable-on-the-first-argument
behaviour, done in-builtin so the bound lists are not split). It rejects
complex-valued `x`, `min`, or `max` via `clip_has_imaginary_part` (emitting
`Clip::ncompl` once through `matsol_warn_once`), handles `Infinity` /
`-Infinity` before any numericalization via `clip_classify_infinity`
(returning `vmax` / `vmin`), then reduces `x`, `min`, `max` to machine doubles
with `clip_to_double_value`. The decision is purely on which side of the bounds
`x` lands: `x < min` returns a copy of `vmin`, `x > max` returns `vmax`,
otherwise the original (symbolic) `x` is returned unchanged. If any of `x`,
`min`, `max` cannot be reduced to a number, the call stays unevaluated so user
DownValues can take over.

**Data structures.** `clip_to_double_value` coerces Integer / Real / BigInt /
MPFR / `Rational[n,d]` directly, and falls back to `numericalize(e,
numeric_machine_spec())` for symbolic constants (so `Clip[Pi]` reduces `Pi` to
~3.14 and yields `1`). Only the comparison is numeric; the returned interior
value preserves the original exact/symbolic `x`.
