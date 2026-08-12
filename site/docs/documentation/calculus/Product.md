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

## Examples (8)

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

**Attributes:** `HoldAll`, `Protected`.

## See also

[Sum](../../calculus/Sum/), [HoldAll](../../expression-information/HoldAll/), [NProduct](../../numerical-calculus/NProduct/), [Pochhammer](../../special-functions/Pochhammer/), [Factorial](../../arithmetic/Factorial/), [Together](../../algebra/Together/), [Factor](../../algebra/Factor/), [QPochhammer](../../special-functions/QPochhammer/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_divisors.c`](https://github.com/stblake/mathilda/blob/main/tests/test_divisors.c)
- Tests: [`tests/test_eigen.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eigen.c)
- Tests: [`tests/test_iter.c`](https://github.com/stblake/mathilda/blob/main/tests/test_iter.c)
