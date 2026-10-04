# Interval

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Interval[{min, max}] represents the range of real values between min and max, inclusive. Interval[{a1,b1}, {a2,b2}, ...] is the union of the ranges. Arithmetic and elementary functions thread through intervals, producing rigorous enclosures; exact endpoints are kept exact.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

Arithmetic threads to a rigorous enclosure

```mathematica
In[1]:= Interval[{1, 2}] + Interval[{3, 4}]
Out[1]= Interval[{4, 6}]
```

The product spans the sign change

```mathematica
In[2]:= Interval[{2, 3}]*Interval[{-1, 1}]
Out[2]= Interval[{-3, 3}]
```

Monotone functions thread; exact endpoints stay exact

```mathematica
In[3]:= Sqrt[Interval[{1, 4}]]
Out[3]= Interval[{1, 2}]
```

The crest at Pi/2 pins the upper endpoint to exact 1

```mathematica
In[4]:= Sin[Interval[{0, Pi}]]
Out[4]= Interval[{0, 1}]
```

Membership test

```mathematica
In[5]:= IntervalMemberQ[Interval[{1, 5}], 3]
Out[5]= True
```

## Implementation notes

**Algorithm.** `Interval[{min, max}, ...]` represents a set of real values — the union of
the given inclusive ranges. `builtin_interval` (`src/interval.c`) canonicalises: it orders
each endpoint pair, sorts the ranges by their lower endpoint, and merges overlapping or
touching ranges, so `Interval[{1,3},{2,5}]` collapses to `Interval[{1,5}]`. It is a plain
`EXPR_FUNCTION` with head `SYM_Interval` (the `Complex`/`Rational` pattern), recognised by
`is_interval`, and marked `Protected` only — deliberately **not** `NumericFunction` or
`Listable`, because an interval is a *set*, not a number threaded over.

Arithmetic and elementary functions thread through intervals to produce **rigorous
enclosures**. The endpoint arithmetic goes through the ordinary evaluator (`Plus`, `Times`,
`Power`, `Sin`, …), so exact endpoints stay exact or symbolic; only an inexact (Real/MPFR)
result is nudged one ULP **outward** by `iv_widen` (`nextafter`, or `mpfr_nextbelow`/
`_nextabove`), which keeps the enclosure sound without any MPFI dependency. The hooks live
at the real dispatch sites: `plus.c`/`times.c`/`power.c` for arithmetic, and
`interval_apply_function` / the general derivative-sign certifier (`interval_thread_call`,
wired into `eval.c` after a builtin returns `NULL`) for monotone and certified-monotone
unary functions. Non-monotone heads such as `Sin`/`Cos` pin an extremum to exact `±1` when a
crest or trough lies inside the range.

**Data structures.** Intervals are ordinary `Expr` trees; endpoint comparison
(`interval_endpoint_cmp`) is done **exactly** — both endpoints are numericalised to 200-bit
MPFR and compared with `mpfr_cmp`, never through the tolerant comparison heads, since a
~1-ULP tolerance there once kept the inner of two candidate endpoints and broke containment.

**Complexity / limits.** Canonicalisation is `O(k log k)` in the number `k` of ranges. The
general certifier is depth-capped (`IV_CERTIFY_MAX_DEPTH`) so oscillatory derivative chains
(Bessel) terminate symbolically rather than looping. `Sqrt` is lowered to `Power[·, 1/2]` at
parse time, so it threads through the power hook, not a `Sqrt` guard. The inclusion property
is fuzzed by `tools/interval_fuzz.py` (`make check-interval`).

**Attributes:** `Protected`.

## References

- R. E. Moore, R. B. Kearfott and M. J. Cloud, *Introduction to Interval Analysis* (SIAM, 2009).
- Source: [`src/interval.c`](https://github.com/stblake/mathilda/blob/main/src/interval.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_elliptic.c`](https://github.com/stblake/mathilda/blob/main/tests/test_elliptic.c)
- Tests: [`tests/test_interval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interval.c)
- Tests: [`tests/test_limit.c`](https://github.com/stblake/mathilda/blob/main/tests/test_limit.c)

## Notes & additional examples

### Notes

`Interval[{min, max}]` represents the inclusive range of reals between `min` and `max`, and
`Interval[{a1,b1}, {a2,b2}, ...]` the union of several ranges; the input is canonicalised
(pairs ordered, ranges sorted and overlapping ones merged). It is a *set*, not a number — it
carries `Protected` only, not `Listable` or `NumericFunction`.

Arithmetic and elementary functions thread through to produce **rigorous enclosures**:
endpoints are computed by the ordinary evaluator, so exact endpoints stay exact (that is why
`Sqrt[Interval[{1,4}]]` is `Interval[{1,2}]` exactly), and only inexact endpoints are nudged
one ULP outward so containment is never violated. Non-monotone functions pin an interior
extremum — `Sin` over `{0, Pi}` reaches its crest at `Pi/2`, giving an upper endpoint of
exactly `1`. Use `IntervalUnion`, `IntervalIntersection` and `IntervalMemberQ` for set
operations.
