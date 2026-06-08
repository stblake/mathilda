# Binomial

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Binomial[n, m]
    gives the binomial coefficient C(n, m) = n! / (m! (n-m)!).
For non-negative integer arguments, computed exactly via GMP's
mpz_bin_uiui. Generalised forms (negative or symbolic n, half-integer
m) reduce through the Gamma functional equation; non-decidable forms
stay unevaluated.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Binomial[10, 3]
Out[1]= 120

In[2]:= Binomial[8.5, -4.2]
Out[2]= 6.04992e-05

In[3]:= Binomial[9/2, 7/2]
Out[3]= 9/2

In[4]:= Binomial[n, 4]
Out[4]= 1/24 n (-3 + n) (-2 + n) (-1 + n)

In[5]:= Binomial[n, n - 1]
Out[5]= n

In[6]:= Binomial[1 + I, 5]
Out[6]= -1/12 - 1/12*I

In[7]:= Binomial[0, 1]
Out[7]= 0
```

## Implementation notes

**Algorithm.** `builtin_binomial` dispatches on argument kind. (1) Integer/integer: it coerces both to `mpz_t` and calls GMP's `mpz_bin_ui` (requires the lower index to fit a `ulong`), returning `0` when the lower index is negative or exceeds a non-negative `n`, and handling negative `n` via the Pascal/upper-negation extension `C(n,k) = (-1)^k C(k-n-1, k)`. (2) Machine reals: the Gamma form `tgamma(n+1)/(tgamma(m+1) tgamma(n-m+1))`. (3) Symmetry/polynomial reduction: if `n - m` evaluates to a small non-negative integer, or `m` itself is a small concrete non-negative integer (`<= 32`), it expands the falling-factorial polynomial via `binomial_polynomial` so that symbolic `n` produces a degree-`m` polynomial that downstream `Expand`/`D` can act on. The `n - m` Subtract is evaluated with arithmetic warnings muted (the exploratory difference may hit a spurious `Power::infy`).

**Data structures.** GMP `mpz_t` for the exact path (results normalised by `expr_bigint_normalize`); symbolic expansions build `Times`/`Plus` trees through `eval_and_free`.

- `Protected`, `Listable`, `NumericFunction`.
- Exact integer/integer path uses GMP (`mpz_bin_ui`), including the

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical Algorithms", on binomial coefficients.
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on the generalized binomial and its polynomial form.
- Source: [`src/numbertheory.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Binomial[50, 25]
Out[1]= 126410606437752
```

```mathematica
In[1]:= Binomial[-1, 3]
Out[1]= -1
```

```mathematica
In[1]:= Binomial[1/2, 3]
Out[1]= 1/16
```

```mathematica
In[1]:= Binomial[n, 2]
Out[1]= 1/2 n (-1 + n)
```

### Notes

Integer arguments give exact coefficients via the falling-factorial product, with
`Binomial[50, 25]` returning the full bigint `126410606437752`. The generalized
definition extends to negative and rational upper arguments: `Binomial[-1, 3] =
-1` and `Binomial[1/2, 3] = 1/16`, using `binomial(x, k) = x(x-1)...(x-k+1)/k!`.
With a symbolic upper argument and a non-negative integer lower argument,
Binomial expands to a polynomial, so `Binomial[n, 2]` becomes `n(n-1)/2`.
