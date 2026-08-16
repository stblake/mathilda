# Binomial

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Binomial[n, m]`**

gives the binomial coefficient C(n, m) = n! / (m! (n-m)!).

**`Gamma[n+1]/(Gamma[m+1] Gamma[n-m+1]): machine precision when an argument`**

<details>
<summary>Notes</summary>

For non-negative integer arguments, computed exactly via GMP's mpz\_bin\_uiui. Generalised forms (negative or symbolic n, half-integer m) reduce through the Gamma functional equation. Rational, real, and complex operands evaluate numerically as is a machine real, and arbitrary precision under N\[Binomial\[..\], p\]; complex arguments reuse the complex Gamma. Non-decidable forms stay symbolic.

</details>

## Examples (21)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (16)

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

In[8]:= Binomial[7/3, 1/5]
Out[8]= Binomial[7/3, 1/5]

In[9]:= N[Binomial[7/3, 1/5]]
Out[9]= 1.33313

In[10]:= N[Binomial[7/3, 1/5], 25]
Out[10]= 1.3331254244650286522359229

In[11]:= Binomial[2.5, 1/5]
Out[11]= 1.34885

In[12]:= Binomial[1/2 + I/3, 1/4]
Out[12]= Binomial[1/2 + 1/3*I, 1/4]

In[13]:= N[Binomial[1/2 + I/3, 1/4]]
Out[13]= 1.08987 + 0.0929283*I

In[14]:= N[Binomial[1/2 + I/3, 1/4], 25]
Out[14]= 1.0898678407199392604353272 + 0.092928304677202434313313055*I

In[15]:= Binomial[1 + I, 5.]
Out[15]= -0.0833333 - 0.0833333*I

In[16]:= Binomial[2. + I, 7 - 3 I]
Out[16]= -75.4683 + 106.815*I
```

### Applications (5)

```mathematica
In[17]:= Binomial[50, 25]
Out[17]= 126410606437752

In[18]:= Binomial[-1, 3]
Out[18]= -1

In[19]:= Binomial[1/2, 3]
Out[19]= 1/16

In[20]:= Binomial[n, 2]
Out[20]= 1/2 n (-1 + n)

In[21]:= Sum[Binomial[4, k] x^k, {k, 0, 4}]
Out[21]= 1 + 4 x + 6 x^2 + 4 x^3 + x^4
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Sum 1/(k(k+1)) to n | 0.283 s | 0.032 s | 1.07 s |
| Sum k^5 to n, closed form | 0.129 s | 0.029 s | 0.116 s |
| Sum 1/k^2 to Infinity | 0.108 s | 0.047 s | 1.38 s |
| Sum r^k to n, symbolic ratio | 0.03 s | 0.034 s | 0.235 s |
| Product (1+1/k) to n | 0.02 s | 0.048 s | 0.29 s |
| Sum Binomial[n,k] over k | 0 s | 0.04 s | 23.8 s |

## Implementation notes

**Algorithm.** `builtin_binomial` dispatches on argument kind. (1) Integer/integer: it coerces both to `mpz_t` and calls GMP's `mpz_bin_ui` (requires the lower index to fit a `ulong`), returning `0` when the lower index is negative or exceeds a non-negative `n`, and handling negative `n` via the Pascal/upper-negation extension `C(n,k) = (-1)^k C(k-n-1, k)`. (2) Machine reals: the Gamma form `tgamma(n+1)/(tgamma(m+1) tgamma(n-m+1))`. (3) Symmetry/polynomial reduction: if `n - m` evaluates to a small non-negative integer, or `m` itself is a small concrete non-negative integer (`<= 32`), it expands the falling-factorial polynomial via `binomial_polynomial` so that symbolic `n` produces a degree-`m` polynomial that downstream `Expand`/`D` can act on. The `n - m` Subtract is evaluated with arithmetic warnings muted (the exploratory difference may hit a spurious `Power::infy`).

**Data structures.** GMP `mpz_t` for the exact path (results normalised by `expr_bigint_normalize`); symbolic expansions build `Times`/`Plus` trees through `eval_and_free`.

- `Protected`, `Listable`, `NumericFunction`.
- Exact integer/integer path uses GMP (`mpz_bin_ui`), including the
  Pascal extension `Binomial[n, m] = (-1)^m Binomial[m-n-1, m]` for
  negative `n`.
- Machine-precision real branch via `tgamma`. Fires when either operand
  is a machine real; the other may be any real-numeric leaf (integer,
  bigint, rational, MPFR) — so `Binomial[2.5, 1/5] → 1.34885`. If an
  operand is symbolic the branch declines and the form stays symbolic
  (`Binomial[2.5, x]`) rather than reading the sibling as `0`.
- Symmetric identity: when `n - m` simplifies to a non-negative
  integer `k ≤ 32`, reduces to `Binomial[n, k]` and expands as a
  falling-factorial polynomial. This catches `Binomial[n, n - 1] → n`,
  `Binomial[9/2, 7/2] → 9/2`, `Binomial[n + 1, n - 1] → n (n + 1)/2`,
  etc.
- Concrete non-negative integer `m ≤ 32` with any other `n` (symbolic,
  rational, complex, …) expands to the falling-factorial polynomial
  `n (n-1) (n-2) ... (n-m+1) / m!`, which the `Times`/`Plus` folders
  then simplify.
- Arbitrary-precision real branch via MPFR (`mpfr_gamma`), evaluating the
  Gamma quotient at the working precision. It leaves the exact rational
  form symbolic — `Binomial[7/3, 1/5]` stays `Binomial[7/3, 1/5]`, as in
  Mathematica — but supplies a value under `N[Binomial[7/3, 1/5], p]`.
  Placed after the integer-`m` branch so an integer `m` keeps the exact
  falling-factorial form. The working precision follows the usual
  contagion rule (minimum precision among the inexact operands, floored
  at machine 53).
- Complex numeric branch: fires when the computation is numeric — an
  inexact operand is present (a machine or MPFR real, or a `Complex[..]`
  with an inexact part) — **and** a complex operand is present, so the
  real-only machine branch has declined. The Gamma quotient is built as an
  expression and evaluated, reusing `Gamma`'s complex kernels (machine
  Lanczos, and the arbitrary-precision Spouge path under MPFR) and the
  complex-arithmetic folders; any exact Gaussian sibling (`1 + I`,
  `7 - 3 I`) is carried along by the `Times`/`Plus` numeric contagion as the
  quotient folds. So `Binomial[1 + I, 5.]`, `Binomial[2. + I, 7 - 3 I]`,
  `N[Binomial[1/2 + I/3, 1/4]]` and `N[Binomial[1/2 + I/3, 1/4], 25]` all
  evaluate. The result is accepted only if it is numeric, so a pole or a
  symbolic operand leaves the form symbolic. A pair of **exact** Gaussians
  with no inexact operand (`Binomial[1 + I, 2 + I]`) stays symbolic.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Times](../../arithmetic/Times/), [Plus](../../arithmetic/Plus/), [Gamma](../../special-functions/Gamma/)

- Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical Algorithms", on binomial coefficients.
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on the generalized binomial and its polynomial form.
- Source: [`src/numbertheory.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_bigint.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bigint.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)

## Notes & additional examples

### Notes

Integer arguments give exact coefficients via the falling-factorial product, with
`Binomial[50, 25]` returning the full bigint `126410606437752`. The generalized
definition extends to negative and rational upper arguments: `Binomial[-1, 3] =
-1` and `Binomial[1/2, 3] = 1/16`, using `binomial(x, k) = x(x-1)...(x-k+1)/k!`.
With a symbolic upper argument and a non-negative integer lower argument,
Binomial expands to a polynomial, so `Binomial[n, 2]` becomes `n(n-1)/2`.
