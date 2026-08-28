# SquareFreeQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SquareFreeQ[expr] gives True if expr is a square-free number or polynomial -- no prime (or irreducible factor) occurs more than once, equivalently MoebiusMu[n] != 0 for an integer -- and False otherwise; the squarefree integers have density 6/Pi^2.`**

**`SquareFreeQ[expr, vars]`**

gives True if expr is square-free with respect to the variables vars.

<details>
<summary>Notes</summary>

Option GaussianIntegers -\> True | False | Automatic switches to Gaussian integers.

</details>

## Examples (18)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (10)

```mathematica
In[1]:= SquareFreeQ[10]
Out[1]= True

In[2]:= SquareFreeQ[4]
Out[2]= False

In[3]:= SquareFreeQ[20]
Out[3]= False

In[4]:= SquareFreeQ[3 + 2 I]
Out[4]= True

In[5]:= SquareFreeQ[2/3]
Out[5]= True

In[6]:= SquareFreeQ[6 + 6 x + x^2]
Out[6]= True

In[7]:= SquareFreeQ[x^3 - x^2 y]
Out[7]= False

In[8]:= SquareFreeQ[x y^2, x]
Out[8]= True

In[9]:= SquareFreeQ[x y^2, y]
Out[9]= False

In[10]:= SquareFreeQ[10^70 + 3]
Out[10]= True
```

### Options (1)

```mathematica
In[11]:= SquareFreeQ[2, GaussianIntegers -> True]
Out[11]= False
```

### Applications (7)

```mathematica
In[12]:= SquareFreeQ[12]
Out[12]= False

In[13]:= SquareFreeQ[30]
Out[13]= True

In[14]:= SquareFreeQ[x^2 - 1]
Out[14]= True

In[15]:= SquareFreeQ[(x - 1)^2 (x + 1)]
Out[15]= False

In[16]:= SquareFreeQ[2]
Out[16]= True

In[17]:= SquareFreeQ[2, GaussianIntegers -> True]
Out[17]= False

In[18]:= SquareFreeQ[x^4 + x^2 + 1]
Out[18]= True
```

## Options & behaviour

Diagnostics:
- `SquareFreeQ[]` emits `SquareFreeQ::argb` and stays unevaluated.
- A non-`Rule` past the optional vars slot (e.g. `SquareFreeQ[1, 2, 3]`)
  emits `SquareFreeQ::nonopt` and stays unevaluated.
- `Modulus -> n` with `n != 0` emits `SquareFreeQ::modnotimpl` and stays
  unevaluated (e.g. `SquareFreeQ[x^2 + 1, Modulus -> 2]`).

## Algorithm

squarefreeq.c -- SquareFreeQ[expr] / SquareFreeQ[expr, vars] / opts.

Always returns True or False on a structurally valid call. Wrong arg count emits `SquareFreeQ::argb` to stderr and returns NULL. Malformed options emit `SquareFreeQ::nonopt` and return NULL.

Algorithms:

```text
  - Integer n: factor with FactorInteger, check every prime exponent <= 1.
              Special-case 0 -> False, +/-1 -> True.
  - Rational p/q: square-free iff both numerator and denominator are.
  - Gaussian integer a + b I (with GaussianIntegers -> True or auto-detected
    from a Complex[Integer, Integer] input): factor N(z) = a^2 + b^2 over Z
    and dispatch by the rational prime's residue mod 4 (see sqfree_gaussian).
  - Polynomial in `vars`: for every var x_i that the polynomial has positive
    degree in, compute PolynomialGCD(p, dp/dx_i); the polynomial is
    square-free iff every such gcd is independent of x_i (degree 0 in x_i).
  - Anything else (Real, symbolic): False -- per Mathematica's "expr is not
    manifestly square free" semantics.
```

The Modulus option is parsed but only Modulus -> 0 is honoured; non-zero values emit `SquareFreeQ::modnotimpl` and return the call unevaluated until a real mod-p polynomial sqfree test is wired in.

## Implementation notes

**Algorithm.** `builtin_squarefreeq` is a `*Q` predicate (always `True`/`False` on a valid call; wrong arg count emits `SquareFreeQ::argb`, malformed options `SquareFreeQ::nonopt`, both returning NULL). `sqfree_dispatch` routes by argument kind: an integer `n` is factored with `FactorInteger` and is square-free iff every prime exponent ≤ 1 (with `0 -> False`, `±1 -> True`); a rational `p/q` iff both numerator and denominator are; a Gaussian integer `a + b I` (with `GaussianIntegers -> True`, or auto-detected from a `Complex[Integer, Integer]`) by factoring the norm `a^2 + b^2` over Z and dispatching on each rational prime's residue mod 4 (`sqfree_gaussian`); a polynomial in `vars` by, for each variable `x_i` of positive degree, computing `PolynomialGCD(p, ∂p/∂x_i)` and requiring it to be degree 0 in `x_i` (independent of the variable). Everything else (Real, symbolic) is `False`.

**Data structures.** `Expr*`; integer/Gaussian work on GMP `mpz_t` via `FactorInteger`; the polynomial path uses the expand + `PolynomialGCD` machinery. `Modulus` is parsed but only `Modulus -> 0` is honoured (non-zero emits `SquareFreeQ::modnotimpl` and leaves the call unevaluated).

- `Protected`. Not `Listable` -- passing a list of inputs treats the list as
  the expression (`SquareFreeQ[{1, 2, 3}]` returns `False`).
- Always returns `True` or `False` on a structurally valid call. For inputs
  that are neither a recognised number nor a manifest polynomial -- reals,
  `Sqrt[2]`, `Sin[x]`, `Pi`, strings -- the result is `False`, never
  symbolic.
- An integer `n` is square-free iff `|n|` has no rational prime factor of
  multiplicity `>= 2`. `0` is not square-free; `+/-1` and `+/-p` (for any
  prime `p`) are.
- A rational `p/q` is square-free iff both `p` and `q` are square-free
  integers.
- A polynomial in `vars` is square-free iff for every variable `x_i` in
  `vars` that the polynomial actually depends on, `PolynomialGCD(p, dp/dx_i)`
  has degree `0` in `x_i`. Implementation routes the derivative through
  the `D` builtin and the gcd through `PolynomialGCD`.
- For `GaussianIntegers -> True` (or `Automatic` on a `Complex[Integer,
  Integer]` input), the test factors `N(z) = a^2 + b^2` over `Z` and
  classifies each rational prime by residue `mod 4`:
  - `p == 2`: the Gaussian prime above `2` is `1 + I`; its multiplicity in
    `z` equals the multiplicity of `2` in `N(z)`.
  - `p ≡ 3 (mod 4)`: `p` itself is a Gaussian prime with `N(p) = p^2`;
    its multiplicity in `z` is half the multiplicity of `p` in `N(z)`.
  - `p ≡ 1 (mod 4)`: `p` splits as `pi * conj(pi)`; ambiguous cases (the
    `e_in_norm == 2` slice) are resolved by computing `pi` via a
    Cornacchia search and stripping the multiplicity directly.
- The `Modulus -> p` option is parsed but only `Modulus -> 0` (the default
  no-modulus path) is wired in; any other value -- non-zero integer,
  non-integer, or symbolic -- emits `SquareFreeQ::modnotimpl` and leaves
  the call unevaluated until a polynomial sqfree-mod-`p` test is added.

**Attributes:** `Protected`.

## References

**See also:** [Symbol](../../expression-information/Symbol/), [List](../../other-advanced/List/), [Complex](../../arithmetic/Complex/), [Pi](../../mathematical-constants/Pi/), [D](../../calculus/D/), [PolynomialGCD](../../algebra/PolynomialGCD/), [Rule](../../assignment-and-rules/Rule/)

- G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — squarefree numbers and their density `6/π²`.
- Source: [`src/poly/squarefreeq.c`](https://github.com/stblake/mathilda/blob/main/src/poly/squarefreeq.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_primenu.c`](https://github.com/stblake/mathilda/blob/main/tests/test_primenu.c)
- Tests: [`tests/test_squarefreeq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_squarefreeq.c)

## Notes & additional examples

### Squarefree numbers

An integer is *squarefree* when no prime divides it more than once — equivalently when
`MoebiusMu[n] ≠ 0`, since `μ(n) = 0` is defined precisely by a repeated prime factor. The
squarefree integers have natural density `6/π² = 1/ζ(2) ≈ 0.6079`. For a polynomial,
squarefreeness means no repeated irreducible factor, detected through `PolynomialGCD` of the
polynomial and its derivative.

### Notes

`SquareFreeQ[expr]` tests an integer or a polynomial for the absence of any
repeated factor. Over the integers it checks the prime factorization; over a
polynomial ring it is decided from `GCD[p, p']` (the polynomial is square-free
exactly when this GCD is constant). The `GaussianIntegers -> True` option moves
the test into `Z[i]`, where rational primes such as `2` can acquire a repeated
factor. A second argument restricts the test to the given variables.
