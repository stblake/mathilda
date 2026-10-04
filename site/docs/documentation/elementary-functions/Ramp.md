# Ramp

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Ramp[x]`**

gives x for x \>= 0 and 0 for x \< 0 -- the positive part of x, and the standard spelling of a rectified linear unit.

<details>
<summary>Notes</summary>

The zero returned for a negative argument carries the argument's own exactness: Ramp\[-1.\] is 0. and Ramp\[-3\] is the exact 0, so a Real list maps to a Real list and an integer list to an integer one. Non-real arguments, and symbolic ones whose sign cannot be decided, are left unevaluated. Ramp is Listable and a NumericFunction.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Ramp[{-1., 0., 2.5}]
Out[1]= {0.0, 0.0, 2.5}

In[2]:= Ramp[{-3, 0, 4}]
Out[2]= {0, 0, 4}

In[3]:= Ramp[{-1/2, 3/4}]
Out[3]= {0, 3/4}

In[4]:= Ramp[1 - Sqrt[2]]
Out[4]= 0

In[5]:= Ramp[1. + 2. I]
Out[5]= Ramp[1.0 + 2.0*I]
```

### Applications (5)

The positive part passes a nonnegative argument through

```mathematica
In[6]:= Ramp[3]
Out[6]= 3
```

A negative argument becomes zero

```mathematica
In[7]:= Ramp[-2]
Out[7]= 0
```

Zero carries the argument's exactness, so a Real stays Real

```mathematica
In[8]:= Ramp[-1.]
Out[8]= 0.0
```

A rectified linear unit applied across a vector

```mathematica
In[9]:= Ramp[{-2, -1., 0, 2.5}]
Out[9]= {0, 0.0, 0, 2.5}
```

An undecidable sign is left unevaluated

```mathematica
In[10]:= Ramp[x]
Out[10]= Ramp[x]
```

## Implementation notes

**Algorithm.** `builtin_ramp` is the positive part max(x, 0), the standard
spelling of a rectified linear unit. It classifies the sign with `ustep_class`:
a non-negative argument is returned unchanged, a negative one becomes zero, and
an undecidable or non-real argument is left unevaluated. The zero returned for a
negative argument carries the **argument's own exactness** — `Ramp[-1.]` is `0.`,
`Ramp[-3]` is the exact `0`, and an MPFR argument returns a zero at its precision —
so a Real vector maps to a Real vector and an integer vector to an integer one
with no mixed-head result (this is why `Ramp` needs no gate on its output where
`Clip` does).

**Data structures.** One `ustep_class` call and at most one `expr_copy`; no
intermediate expressions, unlike `UnitBox`. The ND kernel (`REG_U(Ramp)`) maps a
packed or visible numeric `NDArray` element-wise and preserves the element type,
so the buffer is answered in place.

**Complexity / limits.** `O(1)` per element; `Compile[]` lowers it at scalar and
rank-1 shapes. Like `UnitStep`, `Ramp` threads over an `Interval` argument (it is
non-decreasing, so endpoint threading is a rigorous enclosure). A genuinely
complex argument, or one whose sign cannot be certified, stays symbolic.

- `Listable`, `NumericFunction`, `Protected`.
- The zero returned for a negative argument carries the **argument's own
  exactness**: `Ramp[-1.]` is `0.` and `Ramp[-3]` is the exact `0`. A `Real`
  list therefore maps to a `Real` list and an integer list to an integer one,
  with no mixed-head result -- unlike `Clip`, which returns the *bound* at a
  clipped position and so can put an exact `Integer` into a machine-real answer.
- **Exact symbolic real arguments** are resolved by the same numerical
  certification `UnitStep` uses, so `Ramp[Sqrt[2] - 1]` gives `-1 + Sqrt[2]`
  and `Ramp[1 - Sqrt[2]]` gives `0`.
- Non-real arguments (a `Complex` with non-zero imaginary part) and symbolic
  arguments whose sign cannot be certified are left unevaluated.
- A packed list of `Real`s is handled by a threaded buffer kernel (see
  [`packed-arrays.md`](../packed-arrays/index.md)); an integer buffer materialises, which
  changes speed and not the answer.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Real](../../other-advanced/Real/), [Clip](../../elementary-functions/Clip/), [UnitStep](../../elementary-functions/UnitStep/), [Complex](../../arithmetic/Complex/)

- Source: [`src/piecewise.c`](https://github.com/stblake/mathilda/blob/main/src/piecewise.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_interval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interval.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

`Ramp[x]` is the positive part max(x, 0) — the rectified linear unit (ReLU) of
machine learning, now a single pass where `x UnitStep[x]` once needed two and
produced a mixed Real/Integer product. The zero returned for a negative argument
carries that argument's own exactness, so `Ramp` over a Real vector gives a Real
vector and over an integer vector an integer one, with no mixed-head output.

`Ramp` is non-decreasing, so it threads rigorously over an `Interval` and lowers
under `Compile[]`; its `NDArray` kernel maps a packed buffer element-wise while
preserving the element type.
