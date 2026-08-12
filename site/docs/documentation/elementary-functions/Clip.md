# Clip

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Clip[x]`**

gives x clipped to be between -1 and +1.

**`Clip[x, {min, max}]`**

gives x for min \<= x \<= max, min for x \< min, and max for x \> max.

**`Clip[x, {min, max}, {vmin, vmax}]`**

gives vmin for x \< min and vmax for x \> max.

<details>
<summary>Notes</summary>

Clip threads over lists in its first argument and works at machine or arbitrary precision (via N). Symbolic constants such as Pi are numericalized only to decide which side of the interval x lies on; the original symbolic x is returned unchanged when min \<= x \<= max. Infinity and -Infinity are clipped to the upper and lower replacement values respectively. Clip is not defined for non-real complex values: Clip::ncompl is issued and the call is returned unevaluated. Clip\[a\] for an otherwise undetermined a also stays unevaluated so user-supplied rules can intercept it.

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= Clip[7.5]
Out[1]= 1

In[2]:= Clip[-5/2, {-2, 2}]
Out[2]= -2

In[3]:= Clip[Pi, {-9, 7}, {11, 28}]
Out[3]= Pi

In[4]:= Clip[{-2, 0, 2}]
Out[4]= {-1, 0, 1}

In[5]:= Clip[Infinity]
Out[5]= 1

In[6]:= Clip[Re[2 - 3 I]] + Clip[Im[2 - 3 I]] I
Out[6]= 1 - I

In[7]:= N[Clip[1/11, {1/7, 5}], 50]
Out[7]= 0.142857142857142857142857142857142857142857142857142
```

### Applications (5)

```mathematica
In[1]:= Clip[3.7]
Out[1]= 1
```

```mathematica
In[1]:= Clip[Pi]
Out[1]= 1
```

```mathematica
In[1]:= Clip[{-3, 0.5, 4}, {-1, 1}]
Out[1]= {-1, 0.5, 1}
```

```mathematica
In[1]:= Clip[15, {0, 10}, {-1, 1}]
Out[1]= 1
```

```mathematica
In[1]:= Clip[Infinity, {-2, 2}]
Out[1]= 2
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Clip to [0.25, 0.75] over 4x10^6 | 575 s | 1.95 s | 0.953 s |
| MapThread[Max] over 4x10^6 | 14.8 s | 692 s | 0.772 s |
| MapThread[Min] over 4x10^6 | 14.7 s | 687 s | 0.769 s |
| integer Mod over 4x10^6 | 3.88 s | 0.504 s | 3.28 s |
| a b + a over 4x10^6 | 0.754 s | 1.07 s | 1.41 s |
| a + b over 4x10^6 | 0.383 s | 0.516 s | 0.74 s |

## Implementation notes

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

- `NumericFunction`, `Protected`.
- Threads over a `List` in the first argument: `Clip[{x1, x2, ...}, ...]`
  maps Clip element-wise over the list.  The `{min, max}` and
  `{vmin, vmax}` configuration lists are explicitly **not** threaded
  over -- threading is implemented inside the builtin (not via the
  `Listable` attribute) so the bounds and replacement lists stay intact.
- Symbolic numeric constants (`Pi`, `E`, etc.) are numericalized via
  `N` only to decide which side of the interval `x` lies on; the
  original symbolic `x` is returned when `min <= x <= max`, never
  the numeric approximation.
- `Infinity` and `-Infinity` are handled directly: `Clip[Infinity]`
  yields `vmax` (or the default `1`), `Clip[-Infinity]` yields `vmin`.
- An **infinite bound** means "no limit on that side":
  `Clip[x, {0., Infinity}]` is the positive part and `Clip[x, {-Infinity, 1.}]`
  caps from above only. An infinite bound is never attained, so it cannot put
  its own head into the answer and the result keeps the input's exactness --
  `Clip[{-2., 0.5, 3.}, {0., Infinity}]` is `{0., 0.5, 3.}`, all `Real`. A
  *finite exact* bound beside `Real` data still returns that exact bound where
  it clips (`Clip[{-2., 0.5}, {0, Infinity}]` gives `{0, 0.5}`), which is why
  the two are gated separately on the packed path.
- Complex (non-real) input emits a one-shot `Clip::ncompl` warning and
  the call stays unevaluated.  Use `Re[z]`, `Im[z]` to clip the parts
  separately.
- Symbolic input for which the position cannot be decided
  numerically (e.g. `Clip[a]`) stays unevaluated so user-supplied
  rules can intercept it.

**Attributes:** `NumericFunction`, `Protected`.

## See also

[List](../../other-advanced/List/), [Pi](../../mathematical-constants/Pi/), [E](../../mathematical-constants/E/), [N](../../arithmetic/N/)

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_clip.c`](https://github.com/stblake/mathilda/blob/main/tests/test_clip.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_ndarray_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_reduce.c)

## Notes & additional examples

### Notes

`Clip[x]` saturates `x` to the interval `[-1, 1]`; `Clip[x, {min, max}]`
saturates to `[min, max]`; and `Clip[x, {min, max}, {vmin, vmax}]` returns
`vmin`/`vmax` for out-of-range inputs, giving a piecewise ramp-and-saturate
profile. The first argument threads over lists, so a single call clips a whole
vector. Symbolic constants such as `Pi` are numericalized only to decide which
side of the interval they fall on, and `Infinity`/`-Infinity` clip to the upper
and lower replacement values respectively.
