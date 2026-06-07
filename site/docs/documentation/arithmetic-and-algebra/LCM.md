# LCM

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
LCM[n1, n2, ...]
    gives the least common multiple of the integers ni.
Computed via GMP's mpz_lcm folded across the arguments; sign is
normalised non-negative. Accepts BigInt and Rational inputs.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Flat`, `Listable`, `NumericFunction`, `OneIdentity`, `Orderless`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical Algorithms", on the Euclidean algorithm and least common multiples.
- von zur Gathen & Gerhard, "Modern Computer Algebra", on GCD/LCM relations.
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= LCM[3, 4, 5]
Out[1]= 60
```

```mathematica
In[1]:= LCM[12, 18, 30]
Out[1]= 180
```

```mathematica
In[1]:= LCM[4, 6]
Out[1]= 12
```

```mathematica
In[1]:= LCM[0, 5]
Out[1]= 0
```

### Notes

LCM uses the identity `lcm(a, b) = a*b / gcd(a, b)` and folds across all
arguments, so `LCM[3, 4, 5]` gives `60` and `LCM[12, 18, 30]` gives `180`. The
absorbing convention `LCM[0, n] = 0` holds, matching the fact that zero is the
only common multiple involving zero. Pairwise reduction keeps intermediate values
small, and results promote to GMP bigints when they exceed machine-word range.
