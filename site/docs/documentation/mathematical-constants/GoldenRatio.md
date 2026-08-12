# GoldenRatio

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`GoldenRatio`**

is the golden ratio phi = (1 + Sqrt\[5\])/2, with numerical value ~= 1.61803.

<details>
<summary>Notes</summary>

GoldenRatio is the positive root of x^2 == x + 1. It is a mathematical constant: it has attributes Constant and Protected, NumericQ\[GoldenRatio\] is True, and D\[GoldenRatio, x\] is 0. N\[GoldenRatio, prec\] evaluates it to any precision.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

```mathematica
In[1]:= N[GoldenRatio]
Out[1]= 1.61803
```

```mathematica
In[1]:= N[GoldenRatio, 40]
Out[1]= 1.6180339887498948482045868343656381177203
```

```mathematica
In[1]:= N[GoldenRatio^2 - GoldenRatio - 1, 40]
Out[1]= 0.0
```

```mathematica
In[1]:= FromContinuedFraction[{1, {1}}]
Out[1]= 1/2 (1 + Sqrt[5])
```

```mathematica
In[1]:= Round[N[(GoldenRatio^15 - (1 - GoldenRatio)^15)/Sqrt[5]]]
Out[1]= 610
```

## Implementation notes

- Attributes `Constant`, `Protected`. `Attributes[GoldenRatio] = {Constant,
  Protected}`; the symbol cannot be reassigned.
- Propagated as an exact, unevaluated symbol; `NumericQ[GoldenRatio]` is `True`
  and `D[GoldenRatio, x] = 0`.
- `N[GoldenRatio]` gives the machine value `1.61803`; `N[GoldenRatio, prec]`
  gives any precision, e.g.
  `N[GoldenRatio, 50] = 1.61803398874989484820458683436563811772030917980576`.

**Attributes:** `Constant`, `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)

## Notes & additional examples

### Notes

`GoldenRatio` is `phi = (1 + Sqrt[5])/2`, the positive root of `x^2 == x + 1`;
evaluating that polynomial at `phi` numerically returns `0.0`, confirming the
defining relation. It has the simplest possible continued fraction `[1; 1, 1, ...]`,
so `FromContinuedFraction[{1, {1}}]` recovers it exactly. Through Binet's formula
`Fibonacci[n] == (phi^n - (1 - phi)^n)/Sqrt[5]`, the closed form at `n = 15`
rounds to `610 = Fibonacci[15]`. It is a protected `Constant` (so
`D[GoldenRatio, x]` is `0`) evaluated to any precision by `N`.
