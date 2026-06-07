# RootSum

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
RootSum[Function[t, p[t]], Function[t, body[t]]]
    The formal sum of body[\[Alpha]] over the roots \[Alpha] of
    p[\[Alpha]] == 0.  Held symbolic form, used by the rational
    integrator's NaiveLogPart fallback when the logarithmic part
    cannot be expressed in closed-form real elementary functions.
    Differentiation threads through the body Function:
      D[RootSum[f1, Function[t, body]], x]
        == RootSum[f1, Function[t, D[body, x]]].
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `RootSum[Function[d], Function[body]]` denotes the sum of `body` over the roots of the squarefree polynomial `d`. It is a *held* symbolic head (`ATTR_HOLDALL`); `builtin_rootsum` does not expand the sum but tries one closed-form collapse via `rootsum_try_lagrange`. The recognised identity is the Hermite/Lagrange interpolation collapse

```
  Σ_i  a(α_i) / (d'(α_i) (x − α_i))  ==  a(x) / d(x)
```

valid for squarefree `d` with roots `α_i` and `deg(a) < deg(d)`. This is exactly the form produced by differentiating the log part of a rational integral (`D[RootSum[Function[d], Log[x−#]/d'(#)&], x]`).

The collapse works in `Slot[1]` form: it scans the body for a literal `Plus[x, Times[-1, Slot[1]]]` factor to identify the external variable `x` (`find_x_minus_slot1`), computes `d'` as `D[d, Slot[1]]`, reconstructs `a(#) = body · (x−#) · d'(#)` and simplifies it with `Cancel[Together[...]]` (`simplify_rational`). If the simplified `a` still contains `x` the body did not fit the shape and it bails (`NULL`). Otherwise it substitutes `Slot[1] -> x` in both `a` and `d` (`subst_slot1`) and returns `Together[a(x)/d(x)]`.

**Construction.** `root_make_rootsum` (called from the rational-integration log part, `src/intrat.c`) canonicalises the bound variable into the 1-arg `Function[...Slot[1]...]` form via `substitute_bvar_with_slot`, avoiding leaked context-qualified symbols. D over `RootSum` threads through the body in `src/deriv.c`.

**Limits.** Only the single Lagrange shape is recognised; any other `RootSum` (including numeric expansion over actual roots) is left unevaluated.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Manuel Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed. (Springer, 2005).
- Source: [`src/root.c`](https://github.com/stblake/mathilda/blob/main/src/root.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
