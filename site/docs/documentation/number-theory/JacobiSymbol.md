# JacobiSymbol

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`JacobiSymbol[n, m] gives the Jacobi symbol (n/m); for prime m it is the Legendre symbol, +-1 according to whether n is a quadratic residue modulo m (Euler's criterion) and 0 when m divides n, and it satisfies the law of quadratic reciprocity.`**

<details>
<summary>Notes</summary>

This is the full Kronecker generalisation: the second argument may be even or non-positive and the first may be negative.  Returns -1, 0, or 1.  Listable, and exact via GMP for arbitrary-precision integers.

</details>

## Examples (7)

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

### Applications (2)

```mathematica
In[6]:= JacobiSymbol[2, 7]
Out[6]= 1

In[7]:= JacobiSymbol[3, 7]
Out[7]= -1
```

## Implementation notes

**Algorithm.** `builtin_jacobisymbol` requires two integer-like arguments and returns
`mpz_kronecker(n, m)`, i.e. the full Kronecker-symbol generalisation of the Jacobi symbol:
`m` need not be odd or positive, and `n` may be negative. For prime `m` this is the Legendre
symbol — `+1` if `n` is a non-zero quadratic residue modulo `m`, `-1` for a non-residue, `0`
when `m` divides `n`; for odd composite `m` it is the product of the Legendre symbols over
the prime factors, computed by quadratic reciprocity without factoring `m`. The result is
always `-1`, `0` or `1`.

**Data structures.** Two `mpz_t`, extracted with `expr_to_mpz`; the answer is a machine
`Integer`. There is no ND/packed/`Compile` kernel — `Listable` threading over list arguments
is applied by the evaluator before the builtin is reached.

**Complexity / limits.** `mpz_kronecker` is `O((log m)^2)`, so machine integers and
arbitrary-precision bigints are handled uniformly (e.g. `JacobiSymbol[10^10 + 1,
Prime[1000]]`). Symbolic arguments leave the call unevaluated, and a wrong argument count
emits `JacobiSymbol::argrx`.

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

- K. Ireland and M. Rosen, *A Classical Introduction to Modern Number Theory*, 2nd ed., Springer, 1990 — the Legendre and Jacobi symbols and quadratic reciprocity (Chapter 5).
- G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008.
- K. Ireland and M. Rosen, *A Classical Introduction to Modern Number Theory*, 2nd ed., Springer, 1990 — the Legendre, Jacobi and Kronecker symbols and quadratic reciprocity (Chapter 5).
- Source: [`src/numbertheory/jacobisymbol.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory/jacobisymbol.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_jacobisymbol.c`](https://github.com/stblake/mathilda/blob/main/tests/test_jacobisymbol.c)

## Notes & additional examples

### Quadratic residues and reciprocity

For an odd prime `m`, the Jacobi symbol `(n/m)` is the *Legendre symbol*: `+1` if `n` is a
non-zero quadratic residue modulo `m`, `-1` if it is a non-residue, and `0` if `m ∣ n`. By
*Euler's criterion*, `(n/m) ≡ n^((m-1)/2) (mod m)`. For composite (odd) `m` the Jacobi
symbol is the product of the Legendre symbols over the prime factors, and it obeys the *law
of quadratic reciprocity*, which is what makes it computable in `O(log² n)` steps without
factoring `m` — the same recursion the [`PowerMod`](PowerMod.md) modular square root relies
on. Mathilda returns the full Kronecker generalisation, so `m` may be even or non-positive.
