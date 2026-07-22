# Pochhammer

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Pochhammer[a, n]
    gives the Pochhammer symbol (a)_n = a (a+1) ... (a+n-1) = Gamma[a+n]/Gamma[a].
Exact integer n expands to the product of n linear factors: a polynomial
product for symbolic a, an exact Integer/Rational for numeric a (negative n
gives the reciprocal product). Other numeric arguments evaluate via the
Gamma ratio, reducing exact half-integers to rational multiples of Sqrt[Pi]
and tracking machine or arbitrary (MPFR) precision; machine-precision
complex arguments evaluate too. Other arguments stay symbolic. Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_pochhammer` evaluates the Pochhammer symbol (a)_n = a(a+1)...(a+n-1) = Gamma[a+n]/Gamma[a] by delegating to two existing mechanisms rather than carrying its own numeric kernels. For an exact integer order `n` with `|n| <= POCH_PRODUCT_CAP` (1000) it calls `poch_build_product`, which assembles `Times[a, a+1, ..., a+n-1]` (the reciprocal `Power[Times[...], -1]` for negative `n`) from `Plus[a, k]` factors and runs it through `eval_and_free`; that collapses to an exact Integer/Rational for numeric `a`, preserves MPFR precision and complex arithmetic, and stays a polynomial product for symbolic `a`. Otherwise, when both arguments are `expr_is_numeric_like`, `poch_via_gamma` evaluates `Times[Gamma[a+n], Power[Gamma[a], -1]]` and scans the result with `poch_contains_gamma`, returning `NULL` if a residual `Gamma` head survived (so symbolic inputs stay unevaluated). Short-circuits handle `n == 0` (-> 1), `a == 0` with positive integer `n` (-> 0, before any large factorial), and `Pochhammer[Infinity, n>0]` (-> Infinity). Wrong argument counts emit a `Pochhammer::argrx` diagnostic via `poch_emit_argrx`. Registered in `pochhammer_init` with attributes `Listable`, `NumericFunction`, `Protected`.

- `Pochhammer[a, 0] = 1` for any `a` (including symbolic and `Infinity`).
- Exact integer order `n` (|n| ≤ 1000) expands to the explicit product of

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Abramowitz & Stegun, "Handbook of Mathematical Functions" (1964), §6.1.22 — the Pochhammer symbol (rising factorial).
- NIST Digital Library of Mathematical Functions, §5.2(iii), https://dlmf.nist.gov/5.2 — (a)_n = Gamma(a+n)/Gamma(a).
- Knuth, "The Art of Computer Programming, Vol. 1: Fundamental Algorithms", on rising and falling factorial powers.
- Source: [`src/pochhammer.c`](https://github.com/stblake/mathilda/blob/main/src/pochhammer.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)

## Notes & additional examples

### Worked examples

Exact integer order — a product that collapses to an exact value for numeric
arguments and to a polynomial product for symbolic ones:

```mathematica
In[1]:= Pochhammer[10, 6]
Out[1]= 3603600
```

```mathematica
In[1]:= Pochhammer[n, 5]
Out[1]= n (1 + n) (2 + n) (3 + n) (4 + n)
```

Negative order gives the reciprocal product:

```mathematica
In[1]:= Pochhammer[n, -5]
Out[1]= 1/((-5 + n) (-4 + n) (-3 + n) (-2 + n) (-1 + n))
```

Exact half-integer arguments reduce through the Gamma ratio to rational
multiples of `Sqrt[Pi]`:

```mathematica
In[1]:= Pochhammer[3/2, 1/2]
Out[1]= 2/Sqrt[Pi]
```

Arbitrary precision tracks the input precision (here an exact rational product
numericalised to 50 digits):

```mathematica
In[1]:= N[Pochhammer[1/3, 7], 50]
Out[1]= 505.971650663008687700045724737082761774119798811158
```

Threads element-wise over lists:

```mathematica
In[1]:= Pochhammer[{2, 3, 5, 7, 11}, 3]
Out[1]= {24, 60, 210, 504, 1716}
```

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
