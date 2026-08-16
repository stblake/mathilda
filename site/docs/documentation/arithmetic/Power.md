# Power

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

x ^ y or Power\[x, y\] represents x to the power y. Power is Listable, NumericFunction, and OneIdentity. Integer exponents are reduced exactly (repeated squaring on GMP); Rational and Real exponents evaluate numerically when the base is numeric; Power\[0, 0\] stays Indeterminate; Power\[x, 1/2\] is canonicalised to Sqrt\[x\].

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (8)

```mathematica
In[1]:= Sqrt[45]
Out[1]= 3 Sqrt[5]

In[2]:= (a * b)^2
Out[2]= a^2 b^2

In[3]:= (-1)^(3/2)
Out[3]= -I
```

Floor reduces exponent into [0,1)

```mathematica
In[4]:= (-1)^(7/4) In[4b]:= (-1)^(-1/5) Out[4b]= -(-1)^(4/5)
```

```mathematica
In[5]:= 18^(1/3)
Out[5]= 2^(1/3) 3^(2/3)

In[6]:= 12^(1/3)
Out[6]= 2^(2/3) 3^(1/3)
```

3 and 5 share eff 1/3 -> grouped

```mathematica
In[7]:= 60^(1/3)
Out[7]= 2^(2/3) 15^(1/3)
```

Uniform exps -> stays

```mathematica
In[8]:= 6^(1/3)
Out[8]= 6^(1/3)
```

### Applications (7)

```mathematica
In[9]:= 2^200
Out[9]= 1606938044258990275541962092341162602522202993782792835301376

In[10]:= (1/2)^-5
Out[10]= 32

In[11]:= 27^(2/3)
Out[11]= 9

In[12]:= 0^0
Out[12]= Indeterminate

In[13]:= (3 + 4 I)^10
Out[13]= -9653287 + 1476984*I

In[14]:= Sqrt[-12]
Out[14]= (2*I) Sqrt[3]

In[15]:= N[2^(1/2), 40]
Out[15]= 1.4142135623730950488016887242096980785697
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Fourier 1200000 (mixed radix) | 12.7 s | 9.68 s | 8.53 s |
| ListConvolve 100000 x 2048 | 9.71 s | 1.1 s | 10.2 s |
| ListCorrelate 100000 x 2048 | 9.68 s | 1.1 s | 10.4 s |
| Fourier 262143 (awkward size) | 5.91 s | 5.13 s | 4.42 s |
| Fourier 2^18 (262144) | 4.46 s | 2.7 s | 2.35 s |
| InverseFourier 2^18 | 3.37 s | 2.8 s | 2.39 s |

## Implementation notes

**Algorithm.** `builtin_power` evaluates `Power[base, exp]`. `Power[x]` is x; `Power[b, e1, e2, ...]` is right-associated into `Power[b, Power[e1, e2, ...]]` (right-associative grouping). The two-argument core handles, in order: infinity/`Indeterminate` algebra (`0^Infinity -> 0`, `1^Infinity -> Indeterminate` with message, `Infinity^n` by sign of n, etc.); numeric exact folding (integer/rational/bigint powers via GMP, e.g. exact `2^10`, `(1/2)^3`); inexact Real/MPFR exponentiation; partial radical simplification of `integer^(p/q)` (pulling out perfect-power factors so `Sqrt[8] -> 2 Sqrt[2]`); `(b^m)^n -> b^(m·n)` and product/zero/one identities; and `Sqrt`-style rational-exponent canonicalisation. `Sqrt[x]` is a thin wrapper (`builtin_sqrt`) that rewrites to `Power[x, 1/2]`. Symbolic cases that cannot be reduced return `NULL`, leaving the call unevaluated. `Power` is `ONEIDENTITY | LISTABLE | NUMERICFUNCTION | PROTECTED` (note: not Flat/Orderless — exponentiation is neither associative nor commutative).

**Data structures.** `Expr*` trees; exact integer/bigint exponentiation uses GMP `mpz`, rationals via `make_rational`, and MPFR for high-precision reals. Radical factor extraction works on integer factorisation of the base.

**Complexity / limits.** Integer powers are `O(log exp)` GMP multiplies; radical canonicalisation costs a factorisation of the integer base.

- `Listable`.
- Simplifies integer powers of integers.
- Returns `Overflow[]` if the result exceeds 64-bit integer limits.
- Reduces radicals (e.g., `8^(1/2)` becomes `2*Sqrt[2]`).
- Supports complex results for negative bases (e.g., `(-1)^(1/2)` becomes `I`).
  Higher-power cases for `q == 2` now also reduce: `(-1)^(3/2) → -I`,
  `(-1)^(5/2) → I`, `(-12)^(3/2) → -24 I Sqrt[3]` (the principal-branch
  rule `(-n)^(p/2) = I^p · |n|^(p/2)`).
- For a negative base, the residual `(-1)^(b/q)` exponent is reduced into
  `[0, 1)` by **floor** division of `p/q`, pulling out a `(-1)^a = ±1` sign
  that merges into the coefficient (Mathematica canonical form). This covers
  negative exponents and `|p| ≥ q` alike, e.g. `(-1)^(-1/5) → -(-1)^(4/5)`,
  `(-1)^(-2/3) → -(-1)^(1/3)`, `(-1)^(-7/5) → (-1)^(3/5)`, `(-1)^(5/4) →
  -(-1)^(1/4)`, `(-8)^(-1/3) → -(-1)^(2/3)/2`. Positive bases keep
  truncation toward zero (residual exponent in `(-1, 1)`), unchanged.
  Negative bases other than `-1` with even `q ≥ 4` (e.g. `(-16)^(1/4)`) are
  still left unevaluated.
- Distributes power over product if the exponent is an integer.
- **Nested rational powers compose for any base when `|inner exponent| < 1`.**
  `(B^r)^s → B^(r·s)` holds on the principal branch for *any* complex `B` when
  `r = p/q` is a non-integer rational with `|p| < q` (then `r·Arg(B)` stays in
  `(-π, π]`, so no branch cut is crossed). This works without a positivity
  assumption: `Sqrt[a^(2/3)] → a^(1/3)`, `Sqrt[Sqrt[a]] → a^(1/4)`. Inner
  exponents with `|p| ≥ q` (e.g. `a^(3/2)`) and symbolic inner exponents (e.g.
  `2^a`) still stay unevaluated.
- **Positive numeric coefficient splits out of a mixed `Times` base** under a
  rational power: `(c·w)^(p/q) → c^(p/q) · w^(p/q)` when `c > 0` is a numeric
  rational/integer/real that fully reduces under `q` and `w` is the symbolic
  residual (valid for any `w` since `Arg(c) = 0`). Combined with the nested-power
  rule this gives `Sqrt[(1/27/a)^(2/3)] → 1/3 (1/a)^(1/3)`. Gated to fire only
  when the coefficient genuinely reduces, so `Sqrt[2 Pi]`, `(4 Pi)^(2/3)`,
  `Sqrt[2 Sqrt[3]]` stay nested.
- For `Power[Integer, Rational]` with positive base and positive
  `p/q` exponent, splits the base's prime factorisation into a
  product of distinct-prime powers grouped by reduced effective
  exponent (Mathematica canonical form). Triggers only when the
  resulting form is strictly more informative -- uniform-exponent
  inputs like `6^(1/3)` and `30^(1/3)` keep the compact form.

**Attributes:** `Listable`, `NumericFunction`, `OneIdentity`, `Protected`.

## References

**See also:** [I](../../mathematical-constants/I/), [Times](../../arithmetic/Times/)

- Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical Algorithms", on binary exponentiation.
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on simplification of radical powers.
- Source: [`src/power.c`](https://github.com/stblake/mathilda/blob/main/src/power.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_bigint.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bigint.c)
- Tests: [`tests/test_cherry_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_stress.c)
- Tests: [`tests/test_collect_corpus.c`](https://github.com/stblake/mathilda/blob/main/tests/test_collect_corpus.c)
- Tests: [`tests/test_compile_transforms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_transforms.c)

## Notes & additional examples

### Notes

Integer powers use binary exponentiation and promote to GMP bigints, so `2^200`
is exact. A rational base with a negative integer exponent inverts and raises,
giving `(1/2)^-5 = 32`. Rational exponents trigger perfect-power extraction:
`27^(2/3)` reduces to `9`, while non-extractable cases such as `8^(1/3)` of a
non-cube stay symbolic. The indeterminate form `0^0` evaluates to
`Indeterminate` rather than `1`. Complex bases (Gaussian integers, negative
radicands) are handled in closed form, and irrational powers of numeric bases
evaluate to the requested precision under `N[...]`.
