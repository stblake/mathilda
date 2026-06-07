# Subresultants

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Subresultants[poly1, poly2, var]
    gives the list of principal subresultant coefficients of poly1 and
    poly2 with respect to var.  The list has length
    Min[Exponent[poly1, var], Exponent[poly2, var]] + 1, its first
    element is Resultant[poly1, poly2, var], and the first k entries
    vanish exactly when the polynomials share k roots (with
    multiplicity).  Computed by a subresultant polynomial-remainder
    sequence, or a Sylvester-minor determinant for algebraic
    coefficients.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Subresultants[2x^7 + 3x^3 - 7x + 1, 3x^5 - 17x + 21, x]
Out[1]= {273612691817, 68946901, 1299537, 16641, 0, 9}

In[2]:= Subresultants[(x - 1)(x - 2)(x - 3), (x - 1)(x - 4)(x - 5), x]
Out[2]= {0, 12, -4, 1}

In[3]:= Subresultants[(x - 1)^5 (x - 2)(x - 3), (x - 1)^4 (x - 4)(x - 5), x]
Out[3]= {0, 0, 0, 0, 144, 18, 1}

In[4]:= Subresultants[a x^3 + b x^2 + c x + d, x^3 - 5 b x - 7 a, x]
Out[4]= {-343 a^6 + 245 a^4 b^2 - 49 a^2 b^3 + 147 a^3 b c - 175 a^3 b^2 c + 35 a b^3 c - 70 a^2 b c^2 - 7 a c^3 - 147 a^4 d - 35 a^2 b^2 d + 125 a^2 b^3 d - 25 b^4 d + 21 a b c d + 50 a b^2 c d + 5 b c^2 d - 21 a^2 d^2 - 10 b^2 d^2 - d^3, -7 a^2 b + 25 a^2 b^2 - 5 b^3 + 10 a b c + c^2 - b d, -b, 1}

In[5]:= Length[Subresultants[x^50 + a, x^20 + b, x]]
Out[5]= 21
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
