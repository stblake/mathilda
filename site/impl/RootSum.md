---
references:
  - "Manuel Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed. (Springer, 2005)."
source: src/root.c
---
**Algorithm.** `RootSum[Function[d], Function[body]]` denotes the sum of `body` over the roots of the squarefree polynomial `d`. It is a *held* symbolic head (`ATTR_HOLDALL`); `builtin_rootsum` does not expand the sum but tries one closed-form collapse via `rootsum_try_lagrange`. The recognised identity is the Hermite/Lagrange interpolation collapse

```
  Σ_i  a(α_i) / (d'(α_i) (x − α_i))  ==  a(x) / d(x)
```

valid for squarefree `d` with roots `α_i` and `deg(a) < deg(d)`. This is exactly the form produced by differentiating the log part of a rational integral (`D[RootSum[Function[d], Log[x−#]/d'(#)&], x]`).

The collapse works in `Slot[1]` form: it scans the body for a literal `Plus[x, Times[-1, Slot[1]]]` factor to identify the external variable `x` (`find_x_minus_slot1`), computes `d'` as `D[d, Slot[1]]`, reconstructs `a(#) = body · (x−#) · d'(#)` and simplifies it with `Cancel[Together[...]]` (`simplify_rational`). If the simplified `a` still contains `x` the body did not fit the shape and it bails (`NULL`). Otherwise it substitutes `Slot[1] -> x` in both `a` and `d` (`subst_slot1`) and returns `Together[a(x)/d(x)]`.

**Construction.** `root_make_rootsum` (called from the rational-integration log part, `src/intrat.c`) canonicalises the bound variable into the 1-arg `Function[...Slot[1]...]` form via `substitute_bvar_with_slot`, avoiding leaked context-qualified symbols. D over `RootSum` threads through the body in `src/deriv.c`.

**Limits.** Only the single Lagrange shape is recognised; any other `RootSum` (including numeric expansion over actual roots) is left unevaluated.
