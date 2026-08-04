# DivideBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
DivideBy[x, dx] or x /= dx
    divides x by dx and returns the new value of x.
    x /= dx is equivalent to x = x/dx.

DivideBy has attribute HoldFirst. The first argument x can be a symbol
or a Part expression referring to an existing value. If x has no
assigned value, DivideBy::rvalue is emitted and the expression is left
unevaluated.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= q = 5; q *= 3; q
Out[1]= 15

In[2]:= w = {1., 2., 3.}; w[[3]] /= 4.; w
Out[2]= {1.0, 2.0, 0.75}
```

## Implementation notes

**Attributes:** `HoldFirst`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
