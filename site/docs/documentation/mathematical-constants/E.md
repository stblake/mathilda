# E

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
E
    is the exponential constant e (base of natural logarithms), with
    numerical value ~= 2.71828.
E is a mathematical constant: it has attributes Constant and Protected,
NumericQ[E] is True, and D[E, x] is 0. N[E, prec] evaluates it to any
precision.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Attributes `Constant`, `Protected`. `Attributes[E] = {Constant, Protected}`;

**Attributes:** `Constant`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Log[E^3]
Out[1]= 3
```

```mathematica
In[1]:= N[E, 40]
Out[1]= 2.7182818284590452353602874713526624977572
```

```mathematica
In[1]:= Sum[1/n!, {n, 0, Infinity}]
Out[1]= E
```

```mathematica
In[1]:= Limit[(1 + 1/n)^n, n -> Infinity]
Out[1]= E
```

### Notes

`E` is the exponential constant *e*, the base of the natural logarithm. It is a
protected `Constant`, so `D[E, x]` is `0` and it survives evaluation symbolically
until `N` is applied — `N[E, prec]` returns it to any requested precision via the
MPFR backend. The constant is recognised by the rest of the system, so the
classic limit and series characterisations of *e* both fold back to `E`, and
`Log[E^3]` simplifies to its exponent.
