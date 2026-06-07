# MinimalPolynomial

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
MinimalPolynomial[s, x]
    gives the lowest-degree polynomial in x with integer coefficients,
    positive leading coefficient and content 1, having the algebraic
    number s as a root.  s may be built from rationals, radicals, the
    imaginary unit, roots of unity, and Root[] objects.
MinimalPolynomial[s]
    gives the minimal polynomial as a pure function.
MinimalPolynomial[s, x, Extension -> a]
    gives the characteristic polynomial of s in Q(a) over Q(a).
    Computed by resultant elimination of the radicals; threads over
    lists.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= MinimalPolynomial[Sqrt[2] + Sqrt[3], x]
Out[1]= 1 - 10 x^2 + x^4

In[2]:= MinimalPolynomial[(1 + I)/Sqrt[2], x]
Out[2]= 1 + x^4

In[3]:= MinimalPolynomial[Root[2 #1^3 - 2 #1 + 7 &, 1] + 17, x]
Out[3]= -9785 + 1732 x - 102 x^2 + 2 x^3

In[4]:= MinimalPolynomial[Sqrt[2], x, Extension -> E^(I Pi/4)]
Out[4]= 4 - 4 x^2 + x^4
```

## Implementation notes

- `Listable`, `Protected`. A `List` first argument threads element-wise.
- `s` may be built from integers and rationals, radicals (`Sqrt`,

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
