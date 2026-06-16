# Clip

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Clip[x]
    gives x clipped to be between -1 and +1.
Clip[x, {min, max}]
    gives x for min <= x <= max, min for x < min, and max for x > max.
Clip[x, {min, max}, {vmin, vmax}]
    gives vmin for x < min and vmax for x > max.

Clip threads over lists in its first argument and works at machine
or arbitrary precision (via N). Symbolic constants such as Pi are
numericalized only to decide which side of the interval x lies on;
the original symbolic x is returned unchanged when min <= x <= max.

Infinity and -Infinity are clipped to the upper and lower
replacement values respectively. Clip is not defined for non-real
complex values: Clip::ncompl is issued and the call is returned
unevaluated. Clip[a] for an otherwise undetermined a also stays
unevaluated so user-supplied rules can intercept it.
```

## Examples

All examples below are verified against the current Mathilda build.

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

In[6]:= Clip[2 - 3 I]
Out[6]= Clip[2 - 3*I]

In[7]:= Clip[Re[2 - 3 I]] + Clip[Im[2 - 3 I]] I
Out[7]= 1 - I

In[8]:= N[Clip[1/11, {1/7, 5}], 50]
Out[8]= 0.142857142857142857142857142857142857142857142857142
```

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

**Attributes:** `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)

## Notes & additional examples

### Worked examples

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

### Notes

`Clip[x]` saturates `x` to the interval `[-1, 1]`; `Clip[x, {min, max}]`
saturates to `[min, max]`; and `Clip[x, {min, max}, {vmin, vmax}]` returns
`vmin`/`vmax` for out-of-range inputs, giving a piecewise ramp-and-saturate
profile. The first argument threads over lists, so a single call clips a whole
vector. Symbolic constants such as `Pi` are numericalized only to decide which
side of the interval they fall on, and `Infinity`/`-Infinity` clip to the upper
and lower replacement values respectively.
