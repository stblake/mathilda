# Rationalize

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Rationalize[x]`**

converts an approximate number x to a nearby rational with small denominator.

**`Rationalize[x, dx]`**

yields the rational number with smallest denominator that lies within dx of x.

**`Rationalize[x] yields x unchanged if there is no rational number close enough to x to satisfy |p/q - x| < c/q^2, with c = 10^-4.`**

**`Rationalize[x, dx] works with exact numbers x: the value is first numericalised, then rationalised.`**

**`Rationalize[x, 0] forces conversion of any inexact number x to rational form, using a tolerance derived from the precision of x.`**

<details>
<summary>Notes</summary>

Rationalize threads over compound expressions and Complex\[re, im\], so e.g. Rationalize\[1.2 + 6.7 x\] gives 6/5 + (67 x)/10.

</details>

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (6)

```mathematica
In[1]:= Rationalize[0.5]
Out[1]= 1/2

In[2]:= Rationalize[N[Pi], 10^-4]
Out[2]= 333/106

In[3]:= Rationalize[1.2 + 6.7 x]
Out[3]= 6/5 + 67/10 x
```

Tightening the tolerance recovers the continued-fraction convergents. With the
golden ratio it returns a ratio of consecutive Fibonacci numbers, and with `π`
the celebrated convergent `355/113` and beyond:

```mathematica
In[1]:= Rationalize[N[GoldenRatio], 10^-6]
Out[1]= 1597/987

In[2]:= Rationalize[N[Pi, 40], 10^-10]
Out[2]= 312689/99532
```

`Rationalize[x, 0]` forces conversion of an inexact number using a tolerance
derived from its precision:

```mathematica
In[1]:= Rationalize[N[E], 0]
Out[1]= 325368125/119696244
```

## Algorithm

Mathilda — Rationalize implementation.

Algorithm sketch ---------------- The core of Rationalize is the classical "best rational approximation" problem. Two flavours are needed:

```text
  (A) Rationalize[x, dx] — given a tolerance dx ≥ 0, find the rational
      p/q with the smallest q satisfying |p/q − x| ≤ dx.

      This is solved by the Stern–Brocot / continued-fraction "simplest
      rational in interval" algorithm. We compute the continued fraction
      of x by repeatedly extracting integer parts of the lower and
      upper bound, descending until the integer parts diverge — at which
      point the smallest integer in the lower interval is the simplest
      rational. We use GMP throughout so denominators that exceed
      int64_t (achievable for very small dx, e.g. dx = 0) are handled.

  (B) Rationalize[x] — no tolerance: walk the convergents of the
      continued-fraction expansion of x, returning the first p/q whose
      error satisfies |p/q − x| < 10^-4 / q^2 (the same threshold
      Mathematica uses). If none exist within the precision of x, the
      caller is told (returns false) and is expected to leave x alone.
```

dx = 0 is treated as a thin tolerance derived from the ulp of x — this lets `Rationalize[N[Pi], 0]` produce 245850922/78256779 rather than the bit-exact (and uninteresting) 884279719003555/281474976710656.

Threading --------- Rationalize maps over expression structure. The descent picks one of two strategies per node:

```text
  - Compound expressions that are themselves a NumericQ (e.g. Sqrt[2],
    Pi, Exp[Sqrt[2]]) are first numericalised via the existing N
    pipeline, then rationalised end-to-end. This is what makes
    `Rationalize[Exp[Sqrt[2]], 2^-12]` collapse to 218/53.
  - Otherwise the head and arguments are recursively rationalised
    and then re-evaluated, so e.g. `Rationalize[1.2 + 6.7 x]` becomes
    `6/5 + (67 x)/10` after Plus/Times re-canonicalise.
```

In default mode (no dx) only inexact leaves (EXPR_REAL / EXPR_MPFR) get converted. Exact numbers and symbolic atoms pass through. This matches Mathematica: Rationalize[Pi] = Pi, Rationalize[3.14] = 157/50, ...

Memory ------ The core algorithm allocates / clears its own mpz_t scratch. The caller owns the returned Expr from internal_rationalize_expr / builtin_*. The builtin returns NULL when its input is malformed (caller retains the res Expr); on a successful conversion it returns a fresh Expr and the evaluator releases res itself — the builtin must not free res.

## Implementation notes

**Algorithm.** `builtin_rationalize` dispatches on arity into `internal_rationalize_expr` with one of three modes. `Rationalize[x]` (`RATIONALIZE_DEFAULT`) walks the convergents of the continued-fraction expansion of `x` and returns the first `p/q` whose error satisfies `|p/q − x| < 10^-4 / q^2` (the standard threshold); if none qualify it returns the value unchanged. `Rationalize[x, dx]` with `dx > 0` (`RATIONALIZE_TOLERANCE`) solves the "simplest rational in an interval" problem by the Stern–Brocot / continued-fraction descent: it extracts integer parts of the lower and upper bounds until they diverge, at which point the smallest integer in the lower interval gives the minimal-denominator `p/q` with `|p/q − x| ≤ dx`. `dx == 0` (`RATIONALIZE_ZERO`) is treated as a thin ulp-derived tolerance so e.g. `Rationalize[N[Pi], 0]` yields `245850922/78256779` rather than the bit-exact dyadic. GMP is used throughout so denominators exceeding `int64_t` are handled. A negative or non-numeric `dx` returns `NULL` (left unevaluated).

**Threading.** Rationalize descends through expression structure: a compound subexpression that is itself `NumericQ` (e.g. `Sqrt[2]`, `Exp[Sqrt[2]]`) is numericalised via the `N` pipeline and then rationalised end-to-end; otherwise the head and arguments are recursively rationalised and re-evaluated so `Rationalize[1.2 + 6.7 x]` becomes `6/5 + (67 x)/10`.

**Attributes:** `Protected`.

## References

- Source: [`src/rationalize.c`](https://github.com/stblake/mathilda/blob/main/src/rationalize.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_chop.c`](https://github.com/stblake/mathilda/blob/main/tests/test_chop.c)
- Tests: [`tests/test_rationalize.c`](https://github.com/stblake/mathilda/blob/main/tests/test_rationalize.c)

## Notes & additional examples

### Notes

`Rationalize[x]` finds a nearby rational with small denominator (within `c/q^2`, `c = 10^-4`); the two-argument form `Rationalize[x, dx]` returns the smallest-denominator rational within `dx` of `x`. It threads over compound expressions.
