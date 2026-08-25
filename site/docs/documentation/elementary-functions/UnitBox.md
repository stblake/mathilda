# UnitBox

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`UnitBox[x]`**

gives 1 for -1/2 \<= x \<= 1/2 and 0 otherwise -- the rectangular pulse (box) function, closed at both endpoints.

<details>
<summary>Notes</summary>

The result is always exact. Exact symbolic real arguments are resolved by the same numerical certification UnitStep and Ramp use; non-real or unresolved arguments are left unevaluated. UnitBox is Listable.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= UnitBox[0]
Out[1]= 1

In[2]:= UnitBox[1/2]
Out[2]= 1

In[3]:= UnitBox[-1/2]
Out[3]= 1

In[4]:= UnitBox[0.6]
Out[4]= 0

In[5]:= UnitBox[{-1, -0.5, 0, 0.5, 1}]
Out[5]= {0, 1, 1, 1, 0}

In[6]:= UnitBox[Pi]
Out[6]= 0

In[7]:= UnitBox[x]
Out[7]= UnitBox[x]
```

## Implementation notes

- `Listable`, `NumericFunction`, `Orderless`, `Protected`, matching
  Mathematica. `Orderless` reflects the variadic multidimensional box
  (`UnitBox[x, y, ...]`, like `UnitStep`); this implementation evaluates the
  single-argument pulse.
- The result is **always exact** -- an integer `0` or `1` -- for real numeric
  input, including `Real`/`MPFR` arguments.
- Implemented by reusing `UnitStep`'s sign classifier twice, on `x + 1/2` and
  `1/2 - x`: `x` is in range iff neither shifted value is negative. **Exact
  symbolic real arguments** (`Pi`, `Sqrt[2]`, ...) are therefore resolved by
  the same numerical certification `UnitStep` and `Ramp` use.
- Non-real arguments (a `Complex` with non-zero imaginary part) and
  unresolved symbolic arguments are left unevaluated.
- **Fast paths.** `UnitBox` has a narrowing NDArray kernel (`1` iff
  `-1/2 <= x <= 1/2`, an exact integer, real→int and int→int arms like
  `UnitStep`/`Sign`/`Floor`), so it threads over a visible `NDArray[...]` and
  reads a packed buffer directly; it is on the `AWARE` / `INT64_OK` lists in
  `src/pack.c`, and a packed or int64 array stays packed and exact. It also
  lowers in `Compile[]` (and therefore auto-compiles), scalar and rank-1 array,
  as `(x >= -1/2) (x <= 1/2)` typed as an Integer. The scalar interpreter path
  still reuses `UnitStep`'s sign classifier for exact symbolic-real
  certification.
- Does **not** thread through `Interval` in this version: `Floor`/`Ceiling`
  are the only piecewise functions here that do, because `Interval`
  threading only supports monotone functions, and `UnitBox` (a two-sided box)
  isn't one.

**Attributes:** `Listable`, `NumericFunction`, `Orderless`, `Protected`.

## References

**See also:** [Orderless](../../expression-information/Orderless/), [UnitStep](../../elementary-functions/UnitStep/), [Pi](../../mathematical-constants/Pi/), [Ramp](../../elementary-functions/Ramp/), [Complex](../../arithmetic/Complex/), [Sign](../../arithmetic/Sign/), [Floor](../../arithmetic/Floor/), [Interval](../../other-advanced/Interval/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_unitbox.c`](https://github.com/stblake/mathilda/blob/main/tests/test_unitbox.c)
