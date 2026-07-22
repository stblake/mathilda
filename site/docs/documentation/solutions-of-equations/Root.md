# Root

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Root[Function[t, p[t]], k]
    Represents the k-th root of the univariate polynomial p in the
    variable t. k is canonical: real roots first ascending, then
    complex roots ordered by Re ascending, |Im| ascending, with the
    negative-Im member of each conjugate pair first. N[Root[..]]
    and N[Root[..], prec] return a numerical approximation via a
    companion-matrix + Sturm + Newton pipeline.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= N[Root[Function[#^3 - 2 # - 5], 1], 30]
Out[1]= 2.094551481542326591482386540579

In[2]:= N[Root[Function[#^3 + # + 1], 1], 20]    (* real root first *)
Out[2]= -0.682327803828019327372

In[3]:= N[Root[Function[#^3 + # + 1], 2], 20]    (* conj pair: -Im first *)
Out[3]= 0.341163901914009663686 - 1.16154139999725193609*I

In[4]:= N[Root[Function[#^3 + # + 1], 3], 20]
Out[4]= 0.341163901914009663686 + 1.16154139999725193609*I
```

## Implementation notes

**Algorithm.** `Root` is a held symbolic form. `builtin_root` unconditionally returns `NULL` — the symbol carries `ATTR_HOLDALL | ATTR_PROTECTED`, so by the time the builtin is reached the evaluator has already left the call verbatim, and there is nothing to compute. `Root[Function[t, p[t]], k]` denotes the k-th root of the univariate polynomial p in canonical ordering (real roots first ascending, then complex roots by ascending `Re`, ascending `|Im|`, negative-`Im` member of each conjugate pair first). The useful work happens in callers: `N[Root[..]]` / `N[Root[..], prec]` numericalise via a companion-matrix + Sturm + Newton pipeline (`src/root_numeric.c`); `ToRadicals` (`src/radicals.c`) converts a held `Root` of degree ≤ 4 (or a binomial of higher degree) into closed-form radicals, selecting the k-th root by matching against the numeric value; `D[RootSum[...], x]` threads through the body (`src/deriv.c`); and the rational integrator's `NaiveLogPart` fallback (`src/intrat.c`) constructs `RootSum` nodes when the logarithmic part has no closed-form real expression.

This file also implements the companion head `RootSum`. `builtin_rootsum` is *not* purely held: `rootsum_try_lagrange` recognises the post-differentiation shape `Function[a(#)/(d'(#)(x-#))]` and collapses it via the Hermite/Lagrange interpolation identity `Σ_i a(α_i)/(d'(α_i)(x-α_i)) = a(x)/d(x)` to the closed rational function `a(x)/d(x)`, using `find_x_minus_slot1`, `subst_slot1`, `D`, and `internal_together`/`internal_cancel`. `root_make_rootsum` builds the canonical `Slot[1]` form, rewriting a bound variable to `Slot[1]` throughout.

**Data structures.** `Expr*` trees only. The held form is `Root[Function[poly_in_t], k]`; `RootSum` is `RootSum[Function[poly], Function[body]]` in either `Slot[1]` or 2-arg bound-variable form.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/root.c`](https://github.com/stblake/mathilda/blob/main/src/root.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Root[#^2 - 2 &, 1]
Out[1]= Root[#1^2 - 2 &, 1]

In[2]:= N[Root[#^2 - 2 &, 2], 40]
Out[2]= -1.4142135623730950488016887242096980785697
```

```mathematica
In[1]:= N[Root[#^5 - # - 1 &, 1], 40]
Out[1]= 1.1673039782614186842560458998548421807206
```

```mathematica
In[1]:= {N[Root[#^3 - 2 &, 1], 30], N[Root[#^3 - 2 &, 2], 30], N[Root[#^3 - 2 &, 3], 30]}
Out[1]= {1.259921049894873164767210607278, -0.6299605249474365823836053036392 - 1.09112363597172140356007261419*I, -0.6299605249474365823836053036392 + 1.09112363597172140356007261419*I}
```

```mathematica
In[1]:= {N[Root[#^4 + # + 1 &, 1], 20], N[Root[#^4 + # + 1 &, 2], 20], N[Root[#^4 + # + 1 &, 3], 20], N[Root[#^4 + # + 1 &, 4], 20]}
Out[1]= {-0.72713608449119683998 - 0.430014288329715776416*I, -0.72713608449119683998 + 0.430014288329715776416*I, 0.72713608449119683998 - 0.934099289460529439642*I, 0.72713608449119683998 + 0.934099289460529439642*I}
```

### Notes

`Root[f, k]` is the exact, held representation of the `k`-th root of a
univariate polynomial — including roots like that of `#^5 - # - 1 &`, the
classic example of a quintic with **no radical solution** (Abel–Ruffini), which
`Root` nonetheless names and evaluates to arbitrary precision. The index `k` is
canonical: real roots come first in ascending order, then complex roots ordered
by real part, magnitude of imaginary part, and finally the negative-imaginary
member of each conjugate pair. `N[Root[..], prec]` drives a companion-matrix +
Sturm + Newton pipeline to the requested precision; the `#^3 - 2 &` and
`#^4 + # + 1 &` examples recover the full set of real and complex roots in the
canonical order.
