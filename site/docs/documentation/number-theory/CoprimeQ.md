# CoprimeQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CoprimeQ[n1, n2, ...] yields True if the arguments are pairwise relatively prime -- pairwise GCD equal to 1 -- and False otherwise; two random integers are coprime with probability 6/Pi^2.`**

<details>
<summary>Notes</summary>

Works for machine and BigInt integers.  With GaussianIntegers -\> True, or when any argument is an exact Gaussian integer, coprimality is tested over the Gaussian integers Z\[i\].  Returns False unless the arguments are manifestly coprime; CoprimeQ\[\] is False and CoprimeQ\[n\] is True.  Listable and Orderless.

</details>

## Examples (18)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

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

### Worked examples (10)

```mathematica
In[6]:= CoprimeQ[2^100 - 1, 3^100 - 1]
Out[6]= False

In[7]:= CoprimeQ[2^127 - 1, 2^61 - 1]
Out[7]= True

In[8]:= CoprimeQ[0, 1]
Out[8]= True

In[9]:= CoprimeQ[0, 5]
Out[9]= False

In[10]:= CoprimeQ[6, 35, 143]
Out[10]= True

In[11]:= CoprimeQ[2, 3, 4]
Out[11]= False

In[12]:= CoprimeQ[5 + I, 1 - I]
Out[12]= False

In[13]:= CoprimeQ[2, 5, GaussianIntegers -> True]
Out[13]= True

In[14]:= CoprimeQ[2, 10, GaussianIntegers -> True]
Out[14]= False

In[15]:= CoprimeQ[{1, 2, 3, 4, 5}, 6]
Out[15]= {True, False, False, False, True}
```

### Applications (3)

```mathematica
In[16]:= CoprimeQ[14, 15]
Out[16]= True

In[17]:= CoprimeQ[14, 21]
Out[17]= False

In[18]:= CoprimeQ[6, 35, 143]
Out[18]= True
```

## Implementation notes

**Algorithm.** `builtin_coprimeq` first scans its arguments, lifting an optional
`GaussianIntegers -> True|False` rule out of the numeric operands — because `CoprimeQ`
is `Orderless` the option may sit at any position — and treating any other rule, or a
non-Boolean option value, as malformed (result `False`). With zero operands it returns
`False`, with one operand `True` (there are no pairs). Over the ordinary integers it
extracts every operand with `expr_to_mpz` and tests each unordered pair with `mpz_gcd`,
declaring the set coprime iff every pairwise GCD is `1`. When `GaussianIntegers -> True`
is set, or any operand is an exact non-real Gaussian integer, it switches to `Z[i]`:
`gaussian_pair_coprime` runs the Gaussian Euclidean algorithm with round-to-nearest
division (`coprimeq_round_div` computes `q = floor((2 num + den)/(2 den))`), the remainder
norm strictly decreasing each step, and reports coprime iff the final GCD has norm `1`
(a unit). As a `*Q` predicate it always returns a Boolean, so any operand that is not a
manifest integer or Gaussian integer forces `False`.

**Data structures.** Operands stay as borrowed `Expr*`; each is coerced into a GMP `mpz_t`
(integer path) or a pair of `mpz_t` real/imaginary parts (Gaussian path, via
`coprimeq_to_gaussian`, which accepts a rational integer or a `Complex[a, b]` with
integer-like parts). The Gaussian Euclidean step runs entirely in `mpz_t` scratch
registers. There is no ND/packed/`Compile` kernel: `CoprimeQ` is a Boolean predicate, and
`Listable` threading over list arguments is performed by the evaluator before this builtin
runs.

**Complexity / limits.** The integer path is `O(p^2)` GCDs for `p` operands, each `mpz_gcd`
quasi-linear in the operand size, so machine integers and arbitrary-precision bigints are
handled uniformly. The Gaussian path's Euclidean algorithm always terminates because the
remainder norm decreases. Sign is ignored (GCD uses magnitudes); `CoprimeQ[]` is `False`,
`CoprimeQ[n]` is `True`, and rationals, reals, symbols and malformed options all yield
`False` rather than leaving the call unevaluated.

- Machine integers and GMP bigints, handled uniformly through `mpz_gcd`, so large cases are exact: `CoprimeQ[2^100 - 1, 3^100 - 1]` → `False` (both even), `CoprimeQ[2^127 - 1, 2^61 - 1]` → `True`. Sign is ignored; `GCD(0, n) = |n|`, so `CoprimeQ[0, 1]` → `True` but `CoprimeQ[0, 5]` → `False`.
- More than two arguments are tested *pairwise*: `CoprimeQ[6, 35, 143]` → `True`, while `CoprimeQ[2, 3, 4]` → `False` (2 and 4 share a factor).
- Gaussian integers: with `GaussianIntegers -> True`, or when any argument is an exact Gaussian integer, coprimality is tested over `Z[i]` via the Gaussian Euclidean algorithm (round-to-nearest division). `CoprimeQ[5 + I, 1 - I]` → `False` (both divisible by `1 + I`); `CoprimeQ[2, 5, GaussianIntegers -> True]` → `True`, while `CoprimeQ[2, 10, GaussianIntegers -> True]` → `False`.
- `Orderless`: argument order is irrelevant, and the `GaussianIntegers` option may appear at any position.
- `Listable`: threads element-wise over lists, e.g. `CoprimeQ[{1, 2, 3, 4, 5}, 6]` → `{True, False, False, False, True}`.
- As a `*Q` predicate it always returns a Boolean: `CoprimeQ[]` → `False`, `CoprimeQ[n]` → `True` (no pairs), and anything not a manifestly coprime integer or Gaussian integer — rationals, reals, symbols, malformed options — yields `False` (e.g. `CoprimeQ[a, b]` → `False`).

**Attributes:** `Listable`, `Orderless`, `Protected`.

## References

**See also:** [Orderless](../../expression-information/Orderless/)

- G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — coprimality and the density `6/π²` of coprime pairs.
- Source: [`src/numbertheory/coprimeq.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory/coprimeq.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_coprimeq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_coprimeq.c)
- Tests: [`tests/test_multiplicative_order.c`](https://github.com/stblake/mathilda/blob/main/tests/test_multiplicative_order.c)

## Notes & additional examples

### Relatively prime integers

Two integers are *coprime* (relatively prime) when their greatest common divisor is `1` —
they share no prime factor. `CoprimeQ` tests this **pairwise**, so `CoprimeQ[n1, n2, ...]` is
`True` only when every pair is coprime. A classical density result: two integers drawn at
random are coprime with probability `6/π² = 1/ζ(2) ≈ 0.6079`.
