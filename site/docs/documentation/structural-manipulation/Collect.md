# Collect

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Collect[expr, x]
    expands expr and gathers terms with the same power of x, returning
    a sum of the form Sum[c_k x^k] with each c_k free of x.
Collect[expr, {x1, x2, ...}]
    collects with respect to each xi in turn (nested grouping).
Collect[expr, x, f]
    applies f to each coefficient before re-assembling the sum, useful
    for f = Simplify or f = Factor.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Collect[a x + b y + c x + d y, x]
Out[1]= (a + c) x + b y + d y

In[2]:= Collect[a x^4 + b x^4 + 2 a^2 x - 3 b x + x - 7, x]
Out[2]= -7 + (1 + 2 a^2 - 3 b) x + (a + b) x^4

In[3]:= Collect[a Sqrt[x] + Sqrt[x] + x^(2/3) - c x + 3x - 2b x^(2/3) + 5, x]
Out[3]= 5 + (1 + a) Sqrt[x] + (1 - 2 b) x^(2/3) + (3 - c) x
```

## Implementation notes

- `Protected`.
- Automatically threads over lists, equations, inequalities, and logic functions.
- Effectively writes `expr` as a polynomial in `x` or a fractional power of `x`.
- `Collect[expr, var, h]` applies `h` to the expression that forms the coefficient of each term obtained.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
