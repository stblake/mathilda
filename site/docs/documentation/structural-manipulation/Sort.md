# Sort

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Sort[list] sorts the elements of list into canonical order.
Sort[list, p] sorts using the ordering function p.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Sort[{d, b, c, a}]
Out[1]= {a, b, c, d}

In[2]:= Sort[{Pi, E, 2, 3, 1, Sqrt[2]}, Less]
Out[2]= {1, Sqrt[2], 2, E, 3, Pi}
```

## Implementation notes

- `Protected`.
- Uses an efficient quicksort algorithm.
- Canonical order:
    - Real numbers by numerical value.
    - Complex numbers by real part, then imaginary part magnitude.
    - Strings in dictionary order (lowercase before uppercase).
    - Symbols by name.
    - Expressions by length, then head, then parts depth-first.
- Polynomial order: `x^n` sorts relative to `x`.
- Numeric coefficient stripping: when a `Times` term has a leading numeric factor (Integer, Real, BigInt, MPFR, `Rational[n,d]`, `Complex[re,im]`, or a radical such as `Sqrt[2] = Power[2, 1/2]`), that factor is ignored when computing the term's main factor. As a result, `1 + x^2 + Sqrt[2] x` is canonicalised to `1 + Sqrt[2] x + x^2`, matching Mathematica's main-factor-first ordering.
- Custom ordering function `p` can return `1`, `0`, `-1`, `True`, or `False`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
