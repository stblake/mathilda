# Times

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
x * y * ... or Times[x, y, ...] represents a product of terms.
Times is Flat, Orderless, OneIdentity, Listable, and NumericFunction:
nested Times is auto-flattened, factors are sorted, like factors are
merged into Power, and integer products use exact GMP arithmetic.
Numeric zero collapses the product; a Plus factor is left distributed
(use Expand to distribute).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= 2 * x * 3 * y
Out[1]= 6 x y

In[2]:= I * I
Out[2]= -1

In[3]:= 14/Sqrt[10]
Out[3]= 7 Sqrt[2/5]
```

## Implementation notes

- `Flat`, `Orderless`, `Listable`.
- Combines numeric constants and groups identical bases into `Power` expressions.
- Handles `I` as `Complex[0, 1]`.
- Returns `1` if no arguments are provided.
- Returns `Overflow[]` if integer multiplication overflows or if any argument is `Overflow[]`.
- **Sqrt-coefficient absorption**: a non-trivial rational/integer coefficient

**Attributes:** `Flat`, `Listable`, `NumericFunction`, `OneIdentity`, `Orderless`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical Algorithms", on arbitrary-precision integer multiplication.
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on canonical forms for products.
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= 2^64 * 3
Out[1]= 55340232221128654848
```

```mathematica
In[1]:= (1/3)*(3/7)*7
Out[1]= 1
```

```mathematica
In[1]:= x y x z
Out[1]= x^2 y z
```

### Notes

Times shares the `Flat`, `Orderless`, and `OneIdentity` attributes with Plus, so
factors are flattened and canonically ordered. Repeated symbolic factors are
gathered into powers, turning `x y x z` into `x^2 y z`. Rational factors are
multiplied and reduced to lowest terms, and a product of reciprocals can cancel
exactly to `1`. Integer products promote to GMP bigints automatically when they
exceed machine-word range.
