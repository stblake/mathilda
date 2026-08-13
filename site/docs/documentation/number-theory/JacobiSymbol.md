# JacobiSymbol

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`JacobiSymbol[n, m]`**

gives the Jacobi symbol (n/m).

<details>
<summary>Notes</summary>

For prime m the Jacobi symbol reduces to the Legendre symbol, equal to +-1 according to whether n is a quadratic residue modulo m, and 0 when m divides n.  This is the full Kronecker generalisation: the second argument may be even or non-positive and the first may be negative.  Returns -1, 0, or 1.  Listable, and exact via GMP for arbitrary-precision integers.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= JacobiSymbol[10, 5]
Out[1]= 0

In[2]:= JacobiSymbol[10^10 + 1, Prime[1000]]
Out[2]= 1

In[3]:= JacobiSymbol[7, 6]
Out[3]= 1

In[4]:= JacobiSymbol[{2, 3, 5, 7, 11}, 3]
Out[4]= {-1, 0, -1, 1, -1}

In[5]:= JacobiSymbol[-3, {1, 3, 5, 7}]
Out[5]= {1, 0, -1, 1}
```

## Implementation notes

- `Protected`, `Listable` — threads element-wise over lists and arrays.
- For prime `m` the Jacobi symbol reduces to the Legendre symbol, equal to
  `±1` according to whether `n` is a quadratic residue modulo `m`, and `0`
  when `m` divides `n`.
- Following the Wolfram Language, this is the full Kronecker-symbol
  generalisation: the second argument `m` may be even or non-positive, and
  the first argument `n` may be negative.
- Computed with GMP `mpz_kronecker` in $O((\log m)^2)$ time, so `n` and `m`
  may be arbitrary-precision bignums.
- Non-integer numeric inputs and symbolic arguments flow through
  unevaluated with no diagnostic.
- Diagnostic: `JacobiSymbol::argrx` when called with other than 2 arguments.

**Attributes:** `Listable`, `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_jacobisymbol.c`](https://github.com/stblake/mathilda/blob/main/tests/test_jacobisymbol.c)
