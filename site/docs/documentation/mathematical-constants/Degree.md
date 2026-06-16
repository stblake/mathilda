# Degree

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
Degree
    gives the number of radians in one degree, with numerical value
    Pi/180 (~= 0.0174533).
Multiply by Degree to convert degrees to radians, so 30 Degree is 30
degrees. It is a mathematical constant: it has attributes Constant and
Protected, NumericQ[Degree] is True, and D[Degree, x] is 0. N[Degree,
prec] evaluates it to any precision.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Attributes `Constant`, `Protected`. `Attributes[Degree] = {Constant,

**Attributes:** `Constant`, `Protected`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= N[Degree]
Out[1]= 0.0174533
```

```mathematica
In[1]:= Sin[30 Degree] // N
Out[1]= 0.5
```

```mathematica
In[1]:= N[Tan[60 Degree], 40]
Out[1]= 1.7320508075688772935274463415058723669427
```

```mathematica
In[1]:= N[Degree, 40]
Out[1]= 0.017453292519943295769236907684886127134428
```

### Notes

`Degree` is the mathematical constant `Pi/180`, the number of radians in one degree. Multiply an angle by it to convert degrees to radians, so `30 Degree` is the radian measure of 30 degrees and `Sin[30 Degree]` numericalises to `0.5`. It carries the `Constant` and `Protected` attributes (`D[Degree, x]` is `0`), is recognised by `NumericQ`, and evaluates to arbitrary precision under `N[Degree, prec]` — the last example agrees with `N[Pi/180, 40]` digit for digit.
