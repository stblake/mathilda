# Pochhammer

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Pochhammer[a, n]`**

gives the Pochhammer symbol (a)\_n = a (a+1) ... (a+n-1) = Gamma\[a+n\]/Gamma\[a\].

<details>
<summary>Notes</summary>

Exact integer n expands to the product of n linear factors: a polynomial product for symbolic a, an exact Integer/Rational for numeric a (negative n gives the reciprocal product). Other numeric arguments evaluate via the Gamma ratio, reducing exact half-integers to rational multiples of Sqrt\[Pi\] and tracking machine or arbitrary (MPFR) precision; machine-precision complex arguments evaluate too. Other arguments stay symbolic. Listable.

</details>

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (6)

```mathematica
In[1]:= Pochhammer[10, 6]
Out[1]= 3603600

In[2]:= Pochhammer[n, 5]
Out[2]= n (1 + n) (2 + n) (3 + n) (4 + n)

In[3]:= Pochhammer[n, -5]
Out[3]= 1/((-5 + n) (-4 + n) (-3 + n) (-2 + n) (-1 + n))

In[4]:= Pochhammer[3/2, 1/2]
Out[4]= 2/Sqrt[Pi]

In[5]:= N[Pochhammer[1/3, 7], 50]
Out[5]= 505.971650663008687700045724737082761774119798811158

In[6]:= Pochhammer[{2, 3, 5, 7, 11}, 3]
Out[6]= {24, 60, 210, 504, 1716}
```

## Algorithm

Mathilda -- the Pochhammer symbol (rising factorial).

```text
  Pochhammer[a, n] = a (a+1) ... (a+n-1) = Gamma[a+n] / Gamma[a]
```

The implementation deliberately holds almost no numeric code of its own: two existing, fully-tested mechanisms do the heavy lifting.

```text
  1. Times -- for an exact integer n the symbol expands to the product of
     n linear factors, Times[a, a+1, ..., a+n-1]. Evaluating that product
     collapses to an exact (BigInt / Rational) value for numeric a,
     preserves MPFR precision when a is an EXPR_MPFR, does complex
     arithmetic when a is a Complex[..], and stays a symbolic polynomial
     product for symbolic a. One path covers every kind of a.

  2. Gamma -- for a non-integer (or out-of-range) n the symbol evaluates
     Gamma[a+n]/Gamma[a], reusing the Gamma builtin's exact half-integer
     reductions (-> rational multiples of Sqrt[Pi]), its libm / MPFR real
     paths, and its machine-precision complex Lanczos path -- for free.
```

Mirrors how gamma.c reuses Factorial and how FactorialPower (the falling factorial, numbertheory.c) builds a Times product.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

`builtin_pochhammer` evaluates the Pochhammer symbol (a)_n = a(a+1)...(a+n-1) = Gamma[a+n]/Gamma[a] by delegating to two existing mechanisms rather than carrying its own numeric kernels. For an exact integer order `n` with `|n| <= POCH_PRODUCT_CAP` (1000) it calls `poch_build_product`, which assembles `Times[a, a+1, ..., a+n-1]` (the reciprocal `Power[Times[...], -1]` for negative `n`) from `Plus[a, k]` factors and runs it through `eval_and_free`; that collapses to an exact Integer/Rational for numeric `a`, preserves MPFR precision and complex arithmetic, and stays a polynomial product for symbolic `a`. Otherwise, when both arguments are `expr_is_numeric_like`, `poch_via_gamma` evaluates `Times[Gamma[a+n], Power[Gamma[a], -1]]` and scans the result with `poch_contains_gamma`, returning `NULL` if a residual `Gamma` head survived (so symbolic inputs stay unevaluated). Short-circuits handle `n == 0` (-> 1), `a == 0` with positive integer `n` (-> 0, before any large factorial), and `Pochhammer[Infinity, n>0]` (-> Infinity). Wrong argument counts emit a `Pochhammer::argrx` diagnostic via `poch_emit_argrx`. Registered in `pochhammer_init` with attributes `Listable`, `NumericFunction`, `Protected`.

- `Pochhammer[a, 0] = 1` for any `a` (including symbolic and `Infinity`).
- Exact integer order `n` (|n| ≤ 1000) expands to the explicit product of
  `n` linear factors:
  - Symbolic base → a polynomial product, e.g. `Pochhammer[n, 5] =
    n (1 + n) (2 + n) (3 + n) (4 + n)` and `Pochhammer[x, 4] =
    x (1 + x) (2 + x) (3 + x)`.
  - Numeric base → an exact value, e.g. `Pochhammer[10, 6] = 3603600`,
    `Pochhammer[1, 25] = 25!` (GMP BigInt), `Pochhammer[1/2, 3] = 15/8`.
  - Negative `n` gives the reciprocal product, e.g. `Pochhammer[n, -5] =
    1/((-5 + n) (-4 + n) (-3 + n) (-2 + n) (-1 + n))`,
    `Pochhammer[10, -3] = 1/504`.
  - `Pochhammer[0, n] = 0` for positive integer `n` (short-circuited, so
    even `Pochhammer[0, 1285] = 0`); `Pochhammer[Infinity, n] = Infinity`.
- Other numeric arguments evaluate via the Gamma ratio `Γ(a+n)/Γ(a)`,
  reusing the `Gamma` builtin:
  - Exact half-integers reduce to rational multiples of `Sqrt[Pi]`, e.g.
    `Pochhammer[3/2, 1/2] = 2/Sqrt[Pi]`, `Pochhammer[1/2, 1/2] = 1/Sqrt[Pi]`.
  - Machine-precision real → e.g. `Pochhammer[2.4, 8.5] = 2.31022×10⁶`.
  - Arbitrary precision (MPFR) tracks the input precision, e.g.
    `N[Pochhammer[1/3, 7], 50]` and
    `Pochhammer[1.011111111111000000000000000, 8] = 41552.275849087780380888…`.
  - Machine-precision complex → e.g.
    `Pochhammer[2. + 5 I, 8 I] = 2.13868×10⁻⁶ − 1.42187×10⁻⁵ I`.
- Threads over lists (Listable), e.g.
  `Pochhammer[{2, 3, 5, 7, 11}, 3] = {24, 60, 210, 504, 1716}`.
- All other arguments (e.g. `Pochhammer[a, n]`, `Pochhammer[a, 1/2]`,
  `Pochhammer[1/2, 1/3]`) stay unevaluated. (Derivatives and series, which
  Mathematica expresses through `PolyGamma`, are not yet implemented.)

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Gamma](../../special-functions/Gamma/), [PolyGamma](../../special-functions/PolyGamma/)

- Abramowitz & Stegun, "Handbook of Mathematical Functions" (1964), §6.1.22 — the Pochhammer symbol (rising factorial).
- NIST Digital Library of Mathematical Functions, §5.2(iii), https://dlmf.nist.gov/5.2 — (a)_n = Gamma(a+n)/Gamma(a).
- Knuth, "The Art of Computer Programming, Vol. 1: Fundamental Algorithms", on rising and falling factorial powers.
- Source: [`src/pochhammer.c`](https://github.com/stblake/mathilda/blob/main/src/pochhammer.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_fullsimplify.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fullsimplify.c)

## Notes & additional examples

### Notes

`Pochhammer[a, n]` is the rising factorial (a)ₙ = a (a+1) … (a+n-1) =
Γ(a+n)/Γ(a). The implementation deliberately holds almost no numeric code of
its own. For an exact integer order `n` (with |n| ≤ 1000) it builds the
explicit product `Times[a, a+1, …, a+n-1]` — or its reciprocal for negative
`n` — and evaluates it: that single path yields an exact Integer/Rational for
numeric `a`, preserves arbitrary (MPFR) precision and complex arithmetic, and
stays a symbolic polynomial product for symbolic `a`. For every other numeric
case it evaluates the Gamma ratio `Gamma[a+n]/Gamma[a]`, reusing the `Gamma`
builtin's exact half-integer reductions and its machine, MPFR and complex
numeric paths; a residual-`Gamma` guard keeps genuinely-symbolic inputs
unevaluated.

Useful short-circuits avoid unnecessary work: `Pochhammer[a, 0] = 1` for any
`a`; `Pochhammer[0, n] = 0` for positive integer `n` (so `Pochhammer[0, 1285]`
returns `0` without forming a 1284-term factorial); and
`Pochhammer[Infinity, n] = Infinity` for positive integer `n`. Derivatives and
series — expressible through `PolyGamma` — are not yet
implemented.
