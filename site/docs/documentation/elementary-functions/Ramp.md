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

## Examples (5)

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

## Implementation notes

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

## See also

[Clip](../../elementary-functions/Clip/), [UnitStep](../../elementary-functions/UnitStep/), [Complex](../../arithmetic/Complex/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
