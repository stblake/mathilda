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

This file also implements the companion head `RootSum`. `builtin_rootsum` is *not* purely held: `rootsum_try_lagrange` recognises the post-differentiation shape `Function[a(#)/(d'(#)(x-#))]` and collapses it via the Hermite/Lagrange interpolation identity `Σ_i a(α_i)/(d'(α_i)(x-α_i)) = a(x)/d(x)` to the closed rational function `a(x)/d(x)`, using `find_x_minus_slot1`, `subst_slot1`, `D`, and `internal_together`/`internal_cancel`. `root_make_rootsum` builds the Mathematica-canonical `Slot[1]` form, rewriting a bound variable to `Slot[1]` throughout.

**Data structures.** `Expr*` trees only. The held form is `Root[Function[poly_in_t], k]`; `RootSum` is `RootSum[Function[poly], Function[body]]` in either `Slot[1]` or 2-arg bound-variable form.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/root.c`](https://github.com/stblake/mathilda/blob/main/src/root.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
