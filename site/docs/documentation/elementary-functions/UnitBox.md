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

## Examples (12)

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

### Applications (5)

Inside the unit box, the value is 1

```mathematica
In[8]:= UnitBox[0]
Out[8]= 1
```

The box is closed at both endpoints

```mathematica
In[9]:= UnitBox[1/2]
Out[9]= 1
```

Outside the box, the value is 0

```mathematica
In[10]:= UnitBox[0.7]
Out[10]= 0
```

Listable over a vector of test points

```mathematica
In[11]:= UnitBox[{-1, -1/2, 0, 1/2, 0.7}]
Out[11]= {0, 1, 1, 1, 0}
```

A symbolic argument is left unevaluated

```mathematica
In[12]:= UnitBox[x]
Out[12]= UnitBox[x]
```

## Implementation notes

**Algorithm.** `builtin_unitbox` reuses `ustep_class` twice rather than adding a
second classifier: x is in the box iff neither shifted argument `x + 1/2` nor
`1/2 - x` is certified negative. So `UnitBox[x]` is 1 when both one-sided
`UnitStep`-shaped tests pass, 0 when either shifted argument is negative, and
unevaluated when a side cannot be decided. The box is **closed** at both ends
(`UnitBox[1/2] = 1`), matching `UnitStep[0] = 1`, and the result is always the
exact integer 0 or 1 when determined.

**Data structures.** Each element costs two `Expr` allocations and two
`evaluate()` calls (the shifted arguments `x ± 1/2`), traded for reusing
`ustep_class`'s certification logic instead of duplicating it. The ND kernel
(`ndk_UnitBox_i` double→int64, `ndk_UnitBox_ii` int64→int64, `REG_U`) is
narrowing and lives on `pack.c`'s AWARE + `INT64_OK` list. Unlike `UnitStep`,
`Ramp`, `Round` and `IntegerPart`, `UnitBox` does **not** thread over an
`Interval` argument — the interval machinery encloses only monotone functions and
a two-sided box is not one.

**Complexity / limits.** `O(1)` per element; `Compile[]` lowers the single-argument
box to `CT_INT` at scalar and rank-1 shapes (`UnitBox[0.5]` is `1`, not `1.`). A
non-real or undecidable argument is left unevaluated.

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

**See also:** [Orderless](../../expression-information/Orderless/), [UnitStep](../../elementary-functions/UnitStep/), [Real](../../other-advanced/Real/), [Pi](../../mathematical-constants/Pi/), [Ramp](../../elementary-functions/Ramp/), [Complex](../../arithmetic/Complex/), [Sign](../../arithmetic/Sign/), [Floor](../../arithmetic/Floor/)

- Source: [`src/piecewise.c`](https://github.com/stblake/mathilda/blob/main/src/piecewise.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_unitbox.c`](https://github.com/stblake/mathilda/blob/main/tests/test_unitbox.c)

## Notes & additional examples

### Notes

`UnitBox` is the rectangular pulse: 1 on the closed interval `-1/2 <= x <= 1/2`
and 0 outside it. Both endpoints belong to the box (`UnitBox[1/2] = 1`), matching
the closed-at-zero convention of `UnitStep`, and the result is always the exact
integer 0 or 1 once the argument's position is certified.

Internally the two-sided test reuses `UnitStep`'s one-sided sign certification on
the shifted arguments `x + 1/2` and `1/2 - x`. Because a box is not monotone,
`UnitBox` is the one member of this family that does **not** thread over an
`Interval`; it does carry a narrowing `NDArray` kernel and a `Compile[]` lowering
(`UnitBox[0.5]` compiles to the integer `1`).
