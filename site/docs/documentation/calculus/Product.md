# Product

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Product[f, {i, imax}]`**

gives the product of f for i from 1 to imax.

**`Product[f, {i, imin, imax}], Product[f, {i, imin, imax, di}] and Product[f, {i, {i1, i2, ...}}] use the standard iterator forms; multiple iterators give nested products (an inner bound may depend on an outer index). Product[f, i] gives the indefinite product (anti-quotient). The index is localised (HoldAll). Finite ranges are multiplied out directly; symbolic, indefinite and convergent infinite products are evaluated in exact closed form (n!, Pochhammer, Gamma ratios, base^k, QPochhammer, BarnesG) via a Method polyalgorithm.`**

<details>
<summary>Notes</summary>

Options: Method (Automatic | "Telescoping" | "Rational" | "Geometric" | "QProduct"), VerifyConvergence (default True; a divergent infinite product gives Product::div), GenerateConditions, Assumptions. N\[Product\[...\]\] routes to NProduct.

</details>

## Examples (23)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= Product[k, {k, 1, n}]
Out[1]= Factorial[n]

In[2]:= Product[k + a, {k, 1, n}]
Out[2]= Pochhammer[1 + a, n]

In[3]:= Product[2^k, {k, 1, n}]
Out[3]= 2^(1/2 n (1 + n))

In[4]:= Product[1 + 1/k^2, {k, 1, Infinity}]
Out[4]= Sinh[Pi]/Pi

In[5]:= Product[1 - a q^k, {k, 0, n - 1}]
Out[5]= QPochhammer[a, q, n]

In[6]:= Product[k^k, {k, 1, n}]
Out[6]= Hyperfactorial[n]
```

### Scope (2)

```mathematica
In[7]:= Product[(k^2 - 1)/(k^2 + 1), {k, 2, Infinity}]
Out[7]= Pi Csch[Pi]

In[8]:= Product[(k^3 - 1)/(k^3 + 1), {k, 2, Infinity}]
Out[8]= 2/3
```

### Worked examples (10)

```mathematica
In[9]:= Product[1 - 1/k^2, {k, 2, n}]
Out[9]= (1/2 (1 + n))/n

In[10]:= Product[k, {k, 1, n}]
Out[10]= Factorial[n]

In[11]:= Product[2^(k/2^k), {k, 1, Infinity}]
Out[11]= 4

In[12]:= Product[i^i, {i, 1, n}]
Out[12]= Hyperfactorial[n]

In[13]:= Product[Gamma[i], {i, 1, n-1}]
Out[13]= BarnesG[n]

In[14]:= Product[1 + c/k^2, {k, 1, Infinity}]
Out[14]= Sinh[Pi Sqrt[c]]/(Pi Sqrt[c])

In[15]:= Product[1 + (1/3)^(2^k), {k, 0, Infinity}]
Out[15]= 3/2

In[16]:= Product[Cos[Pi/2^(k+1)], {k, 1, Infinity}]
Out[16]= 2/Pi

In[17]:= Product[Cos[x/2^k], {k, 1, Infinity}]
Out[17]= Sin[x]/x

In[18]:= Product[1/(1 - Prime[i]^-s)]
Out[18]= Product[1/(1 - Prime[i]^(-s))]
```

### Applications (5)

A symbolic finite product is a factorial

```mathematica
In[19]:= Product[k, {k, 1, n}]
Out[19]= Factorial[n]
```

A finite numeric range is multiplied out

```mathematica
In[20]:= Product[k^2, {k, 1, 5}]
Out[20]= 14400
```

A telescoping rational product

```mathematica
In[21]:= Product[(k + 1)/k, {k, 1, n}]
Out[21]= 1 + n
```

A polynomial-exponential (geometric) product

```mathematica
In[22]:= Product[2^k, {k, 1, n}]
Out[22]= 2^(1/2 n (1 + n))
```

A convergent infinite product

```mathematica
In[23]:= Product[1 - 1/k^2, {k, 2, Infinity}]
Out[23]= 1/2
```

## Algorithm

product.c -- Product dispatcher for Mathilda.

```text
The multiplicative analogue of Sum (src/sum/sum.c).  Product is HoldAll: the
```

product variable and bounds must be held so that the iterator is not prematurely evaluated against an outer binding (exactly as Sum/Table/Do hold their iterator specs).

Responsibilities of this file (Stage 0):

```text
  - strip trailing options (Method -> "...", VerifyConvergence -> ..., etc.);
  - rewrite multiple iterators Product[f, s1, ..., sk] into nested single-spec
    products (outer-depends-on-inner bounds come for free);
  - finite explicit expansion: when a range resolves to a finite span of
    integers, or the spec iterates an explicit list, bind the variable and
    fold the evaluated terms with Times (an empty product is 1);
  - otherwise (symbolic bounds, Infinity, or the indefinite form Product[f,i])
    run a Method cascade over the context-qualified sub-algorithms
    Product`Telescoping, Product`Rational, Product`Geometric, Product`QProduct.
    Each sub-builtin returns the closed form (definite:
    Product`M[f,i,imin,imax]; indefinite: Product`M[f,i]) or comes back
    unevaluated to signal "fall through".  When all stages fall through the
    Product[...] is returned unevaluated (held).
```

Adding a later stage is purely additive: a new src/product/product_*.c file, one try_* line in the cascade, and one *_init() call in product_init().

Memory contract: builtin_product takes ownership of res but must not free it

```text
(the evaluator owns it).  Every Expr* allocated here is freed on all paths.
```

## Implementation notes

**Algorithm.** `builtin_product` is the multiplicative analogue of `Sum`, and is
`HoldAll` so the iterator is not evaluated against an outer binding. It strips
trailing options (`Method`, `VerifyConvergence`), rewrites multiple iterators
`Product[f, s1, ..., sk]` into nested single-spec products, then dispatches on
the one spec in `product_one_spec`. A finite numeric range or an explicit list
is multiplied out directly: `expand_list`/`expand_range` bind the index and fold
the evaluated terms with `Times` (an empty product is `1`). Before the slow
evaluate-per-factor loop, a unit-step integer range whose body is inexact at the
first index takes a compiled machine multiply–accumulate (`product_try_compiled`,
via auto-compilation); an exact body (e.g. `n!`) stays on the interpreter so it
keeps exact bignum arithmetic.

**Data structures.** `Expr` trees, the `IterSpec` lattice parser
(`iter_spec_parse_lattice`), and the iterator shadow/restore pair
(`iter_spec_shadow`/`iter_spec_restore`) that localises the bound index. For
symbolic bounds, `Infinity`, or the indefinite form `Product[f, i]`, it runs a
`Method` cascade over context-qualified sub-builtins
(`Product\`Telescoping`, `Product\`Rational`, `Product\`Geometric`,
`Product\`QProduct`, then `Product\`Special`, `Product\`Cantor`,
`Product\`Viete`, `Product\`EulerPrime`, `Product\`RationalInfinite`,
`Product\`BesselZero`, `Product\`Infinite`, `Product\`LogSum`), each returning a
closed form or coming back unevaluated to signal "fall through"
(`result_is_unresolved`). The cascade is ordered cheapest/most-specific first so
the nicest closed form wins.

**Complexity / limits.** A closed-form stage is independent of the span width; a
finite enumeration is linear in the number of factors (guarded by
`PRODUCT_MAX_FINITE_TERMS`, 10^8). The closed-form stages assume a unit step, so
a non-unit step (`{i, 1, n, 2}`) with symbolic/over-wide bounds is left held
rather than given a wrong step-1 form. A finite range whose body hides an
index-dependent predicate (`EvenQ[k]`, `PrimeQ[k]`, ...) is forced to enumerate,
since symbolic evaluation would collapse the predicate and telescope the wrong
factor. When every stage falls through, `Product[...]` is returned unevaluated.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Sum](../../calculus/Sum/), [HoldAll](../../expression-information/HoldAll/), [NProduct](../../numerical-calculus/NProduct/), [Pochhammer](../../special-functions/Pochhammer/), [Factorial](../../arithmetic/Factorial/), [Together](../../algebra/Together/), [Factor](../../algebra/Factor/), [QPochhammer](../../special-functions/QPochhammer/)

- Source: [`src/product/product.c`](https://github.com/stblake/mathilda/blob/main/src/product/product.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_divisors.c`](https://github.com/stblake/mathilda/blob/main/tests/test_divisors.c)
- Tests: [`tests/test_eigen.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eigen.c)

## Notes & additional examples

### Notes

`Product` is the multiplicative analogue of `Sum`. It is `HoldAll`, so the index
is localised and the iterator bounds are not evaluated against an outer binding.
A finite numeric range (or an explicit list of values) is multiplied out
directly, with an empty product giving `1`; a symbolic, indefinite, or convergent
infinite product is handed to a closed-form method cascade.

The cascade tries, cheapest-first, the telescoping (Gamma-free rational),
rational (Pochhammer / Gamma), geometric (`base^k`) and q-product families, plus
several infinite-product specialists. The method can be pinned with
`Method -> "Telescoping" | "Rational" | "Geometric" | "QProduct"`, and
convergence testing for infinite products can be disabled with
`VerifyConvergence -> False`. Multiple iterators `Product[f, s1, s2]` form nested
products, so an inner bound may depend on an outer index. When no method applies
the `Product[...]` is returned unevaluated.
