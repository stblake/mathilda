# Expand

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Expand[expr] expands out products and powers in expr.`**

**`Expand[expr, patt] leaves unexpanded any parts of expr that are free of the pattern patt.`**

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= Expand[(x+3)(x+2)]
Out[1]= 6 + 5 x + x^2

In[2]:= Expand[(x+y)^2 (x-y)^2]
Out[2]= x^4 - 2 x^2 y^2 + y^4

In[3]:= Expand[(x+1)^2 + (y+1)^2, x]
Out[3]= 1 + 2 x + x^2 + (1 + y)^2

In[4]:= Expand[(a+b)(x[1]+x[2])^2, x[_]]
Out[4]= 2 (a + b) x[1] x[2] + (a + b) x[1]^2 + (a + b) x[2]^2

In[5]:= Length[Expand[(1 + x + 5 x^3 + 8 x^17)^341]]
Out[5]= 5758

In[6]:= Expand[(a+b+c+d+e+f+g)^60]
Out[6]= Overflow[]
```

### Applications (6)

A binomial power: the coefficients are the binomial coefficients

```mathematica
In[7]:= Expand[(x + 1)^3]
Out[7]= 1 + 3 x + 3 x^2 + x^3
```

A product of sums distributes into every pairwise term

```mathematica
In[8]:= Expand[(a + b)(c + d)]
Out[8]= a c + b c + a d + b d
```

Several variables: the multinomial, not just the binomial, case

```mathematica
In[9]:= Expand[(1 + x + y)^3]
Out[9]= 1 + 3 x + 3 x^2 + x^3 + 3 y + 6 x y + 3 x^2 y + 3 y^2 + 3 x y^2 + y^3
```

Numeric coefficients are folded, so terms collect and cancel

```mathematica
In[10]:= Expand[(x + 2)^2 (x - 1)]
Out[10]= -4 + 3 x^2 + x^3
```

The work is polynomial in the exponent, not exponential

```mathematica
In[11]:= Expand[(1 + x)^10]
Out[11]= 1 + 10 x + 45 x^2 + 120 x^3 + 210 x^4 + 252 x^5 + 210 x^6 + 120 x^7 + 45 x^8 + 10 x^9 + x^10
```

Only the numerator expands: Expand leaves denominators alone

```mathematica
In[12]:= Expand[(x + 1)^2/(y + 1)]
Out[12]= (1 + 2 x + x^2)/(1 + y)
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Factor x^120 - 1 | 7.01 s | 0.044 s | 3.85 s |
| Discriminant of deg 20 | 2.51 s | 0.068 s | 0.182 s |
| Factor dense degree 60 | 2.2 s | 0.599 s | 26.5 s |
| Factor product of 16 quadratics | 1.29 s | 0.722 s | 27.6 s |
| Factor sparse degree 60 | 0.532 s | 0.516 s | 9.48 s |
| Factor product of 8 quadratics | 0.46 s | 0.191 s | 3.31 s |

## Implementation notes

**Algorithm.** `builtin_expand` (in `src/expand.c`) calls `expr_expand_patt(e, patt)`, a structural recursion that distributes products over sums. The dispatch by head:

- **Plus** — expand each summand and rebuild via `eval_and_free`.
- **Times** — expand each factor, then multiply them pairwise with `multiply_all`, a divide-and-conquer fold whose leaf operation `multiply_two` handles the four `Plus×Plus / Plus×atom / atom×Plus / atom×atom` cases, producing every cross term (`a_i · b_j`) and re-summing them.
- **Power[base, k]** with `k` a positive integer `< 100` — expand the base, then `power_expand` raises it by **repeated squaring** (binary exponentiation), each multiply going through the distributing `multiply_two`. (Note: there is no explicit binomial-coefficient formula; the binomial expansion of `(a+b)^k` falls out of the iterated distribution.) Negative or non-integer exponents are left untouched.
- **List / equations / inequalities / And / Or / Not** — threads into each argument (passing operator-symbol slots of `Inequality` through verbatim).

A two-argument `Expand[expr, patt]` only expands subexpressions that *contain* `patt` (gated by `expr_contains_patt`, which uses the pattern matcher); everything else is copied unchanged. Expand does **not** descend into the arguments of arbitrary function heads.

**Data structures.** Plain `Expr*` trees; transient `Expr**` argument buffers per node. All arithmetic rebuilding goes back through the evaluator (`eval_and_free`) so `Plus`/`Times` canonicalisation (Flat/Orderless, like-term collection) applies after each step.

**Complexity / limits.** Worst case is the full cross-product blow-up of distribution: expanding `(x_1+...+x_m)^k` produces `O(m^k)`-sized output, and `Power` is capped at exponent `< 100` to bound it.

- `Protected`.
- Works only on positive integer powers and distributes products over sums.
- Threads over equations, inequalities, and lists.
- **Polynomials over the rationals are expanded through FLINT** (packed
  `fmpq_mpoly` arithmetic), which distributes and collects like terms orders of
  magnitude faster than the generic tree multiplier on dense, high-degree
  inputs — e.g. `(1 + x + 5 x^3 + 8 x^17)^341` (degree 5797, 5758 terms) expands
  in a few hundredths of a second. Non-polynomial parts (transcendental heads,
  symbolic/fractional exponents, inexact coefficients) fall back to the generic
  binary-splitting distributor and repeated-squaring power expander.
- In a **mixed product** — one where a non-polynomial factor (e.g. `Log[x]`)
  would otherwise force the whole product onto the generic path — the
  polynomial-over-`Q` factors are still multiplied together through FLINT, and
  only the non-polynomial factors are distributed generically. So
  `Expand[Log[x] (1+y)^300 (1-y)^300]` collapses the two degree-300 factors in
  packed arithmetic rather than distributing 90 601 terms by hand.
- `ExpandNumerator` and `ExpandDenominator` inherit the same acceleration (they
  expand their numerator/denominator polynomial through the same path).
- **Never silently declines.** `Expand` expands anything that fits in memory; an
  expansion whose (estimated) size exceeds the memory ceiling returns
  `Overflow[]` rather than leaving the expression unevaluated. The size estimate
  is a Newton-box upper bound — the tighter of the multinomial term count and
  the per-variable degree box — so univariate and low-dimensional expansions of
  *any* degree are recognised as cheap and always run; only a genuinely
  high-dimensional combinatorial blow-up (e.g. `(a+b+c+d+e+f+g)^60`) overflows.
- `Expand[expr, patt]` leaves unexpanded any parts of `expr` that are free of the pattern `patt`. Inside a product, pattern-free factors are carried along as an unexpanded coefficient rather than being distributed.

**Attributes:** `Protected`.

## References

**See also:** [ExpandNumerator](../../algebra/ExpandNumerator/), [ExpandDenominator](../../algebra/ExpandDenominator/)

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), Ch. 3 (normal forms and the distributive expansion of polynomials).
- Source: [`src/expand.c`](https://github.com/stblake/mathilda/blob/main/src/expand.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_bernoullib.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bernoullib.c)
- Tests: [`tests/test_characteristicpolynomial.c`](https://github.com/stblake/mathilda/blob/main/tests/test_characteristicpolynomial.c)
- Tests: [`tests/test_collect_corpus.c`](https://github.com/stblake/mathilda/blob/main/tests/test_collect_corpus.c)
- Tests: [`tests/test_deriv.c`](https://github.com/stblake/mathilda/blob/main/tests/test_deriv.c)

## Notes & additional examples

### Notes

`Expand` applies the distributive law to products and integer powers,
producing a flat sum of monomials in canonical order (ascending total
degree in the leading variable). Like terms are combined automatically, so
`(x + 2)^2 (x - 1)` collapses the `x^1` coefficient to zero and it drops out
of the result. `Expand` only multiplies out — it does not factor or cancel —
and a second argument `Expand[expr, patt]` leaves alone any parts free of the
pattern `patt`.
