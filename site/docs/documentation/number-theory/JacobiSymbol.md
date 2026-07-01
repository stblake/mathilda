# JacobiSymbol

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
JacobiSymbol[n, m]
    gives the Jacobi symbol (n/m).

For prime m the Jacobi symbol reduces to the Legendre symbol, equal to +-1 according to whether n is a quadratic residue modulo m, and 0 when m divides n.  This is the full Kronecker generalisation: the second argument may be even or non-positive and the first may be negative.  Returns -1, 0, or 1.  Listable, and exact via GMP for arbitrary-precision integers.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= JacobiSymbol[10, 5]
Out[1]= 0

In[2]:= Table[JacobiSymbol[n, m], {n, 0, 10}, {m, 1, n, 2}]
Out[2]= {{}, {1}, {1}, {1, 0}, {1, 1}, {1, -1, 0}, {1, 0, 1}, {1, 1, -1, 0}, {1, -1, -1, 1}, {1, 0, 1, 1, 0}, {1, 1, 0, -1, 1}}

In[3]:= JacobiSymbol[10^10 + 1, Prime[1000]]
Out[3]= 1

In[4]:= JacobiSymbol[7, 6]
Out[4]= 1

In[5]:= JacobiSymbol[{2, 3, 5, 7, 11}, 3]
Out[5]= {-1, 0, -1, 1, -1}

In[6]:= JacobiSymbol[-3, {1, 3, 5, 7}]
Out[6]= {1, 0, -1, 1}
```

## Implementation notes

- `Protected`, `Listable` — threads element-wise over lists and arrays.
- For prime `m` the Jacobi symbol reduces to the Legendre symbol, equal to

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
