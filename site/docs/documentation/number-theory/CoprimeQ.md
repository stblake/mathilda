# CoprimeQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
CoprimeQ[n1, n2, ...]
    yields True if the arguments are pairwise relatively prime, and
    False otherwise.
Integers are relatively prime when their GCD is 1.  Works for machine
and BigInt integers.  With GaussianIntegers -> True, or when any
argument is an exact Gaussian integer, coprimality is tested over the
Gaussian integers Z[i].  Returns False unless the arguments are
manifestly coprime; CoprimeQ[] is False and CoprimeQ[n] is True.
Listable and Orderless.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= CoprimeQ[8, 11]
Out[1]= True

In[2]:= CoprimeQ[2, 4]
Out[2]= False

In[3]:= CoprimeQ[2, 3, -5, 7]
Out[3]= True

In[4]:= CoprimeQ[5 + I, 1 - I]
Out[4]= False

In[5]:= CoprimeQ[{1, 2, 3, 4, 5}, 6]
Out[5]= {True, False, False, False, True}
```

## Implementation notes

- Machine integers and GMP bigints, handled uniformly through `mpz_gcd`, so large cases are exact: `CoprimeQ[2^100 - 1, 3^100 - 1]` → `False` (both even), `CoprimeQ[2^127 - 1, 2^61 - 1]` → `True`. Sign is ignored; `GCD(0, n) = |n|`, so `CoprimeQ[0, 1]` → `True` but `CoprimeQ[0, 5]` → `False`.
- More than two arguments are tested *pairwise*: `CoprimeQ[6, 35, 143]` → `True`, while `CoprimeQ[2, 3, 4]` → `False` (2 and 4 share a factor).
- Gaussian integers: with `GaussianIntegers -> True`, or when any argument is an exact Gaussian integer, coprimality is tested over `Z[i]` via the Gaussian Euclidean algorithm (round-to-nearest division). `CoprimeQ[5 + I, 1 - I]` → `False` (both divisible by `1 + I`); `CoprimeQ[2, 5, GaussianIntegers -> True]` → `True`, while `CoprimeQ[2, 10, GaussianIntegers -> True]` → `False`.
- `Orderless`: argument order is irrelevant, and the `GaussianIntegers` option may appear at any position.
- `Listable`: threads element-wise over lists, e.g. `CoprimeQ[{1, 2, 3, 4, 5}, 6]` → `{True, False, False, False, True}`.
- As a `*Q` predicate it always returns a Boolean: `CoprimeQ[]` → `False`, `CoprimeQ[n]` → `True` (no pairs), and anything not a manifestly coprime integer or Gaussian integer — rationals, reals, symbols, malformed options — yields `False` (e.g. `CoprimeQ[a, b]` → `False`).

**Attributes:** `Listable`, `Orderless`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
