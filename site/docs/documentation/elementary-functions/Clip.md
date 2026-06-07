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

- `NumericFunction`, `Protected`.
- Threads over a `List` in the first argument: `Clip[{x1, x2, ...}, ...]`

**Attributes:** `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
