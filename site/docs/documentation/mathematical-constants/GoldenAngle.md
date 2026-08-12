# GoldenAngle

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GoldenAngle`**

is the golden angle (3 - Sqrt\[5\]) Pi = 2 Pi / GoldenRatio^2, with numerical value ~= 2.39996 radians (~= 137.5 degrees).

**`N[GoldenAngle, prec] evaluates it to any precision.`**

<details>
<summary>Notes</summary>

GoldenAngle is a mathematical constant: it has attributes Constant and Protected, NumericQ\[GoldenAngle\] is True, and D\[GoldenAngle, x\] is 0.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= N[GoldenAngle]
Out[1]= 2.39996
```

```mathematica
In[1]:= N[GoldenAngle, 40]
Out[1]= 2.399963229728653322231555506633613853125
```

```mathematica
In[1]:= N[GoldenAngle/Degree, 30]
Out[1]= 137.5077640500378546463487396284
```

```mathematica
In[1]:= N[GoldenAngle - (3 - Sqrt[5]) Pi, 40]
Out[1]= 0.0
```

## Implementation notes

- Attributes `Constant`, `Protected`. `Attributes[GoldenAngle] = {Constant,
  Protected}`; the symbol cannot be reassigned.
- Propagated as an exact, unevaluated symbol; `NumericQ[GoldenAngle]` is `True`
  and `D[GoldenAngle, x] = 0`.
- `N[GoldenAngle]` gives the machine value `2.39996`; `N[GoldenAngle, prec]`
  gives any precision (computed in MPFR from the closed form `(3 - Sqrt[5])
  Pi`), e.g.
  `N[GoldenAngle, 50] = 2.3999632297286533222315555066336138531249990110581`.

**Attributes:** `Constant`, `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)

## Notes & additional examples

### Notes

`GoldenAngle` is `(3 - Sqrt[5]) Pi = 2 Pi / GoldenRatio^2`, the angle that
divides a full turn in the golden ratio — about `137.5` degrees, the divergence
angle that governs optimal phyllotactic spiral packing in plants. Dividing by
`Degree` shows the familiar `137.5...`, and subtracting the closed form
`(3 - Sqrt[5]) Pi` numerically returns `0.0`. It is a protected `Constant`
(so `D[GoldenAngle, x]` is `0`) that `N` evaluates to any precision.
