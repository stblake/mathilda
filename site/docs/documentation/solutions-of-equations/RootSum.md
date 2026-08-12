# RootSum

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RootSum[Function[t, p[t]], Function[t, body[t]]]`**

The formal sum of body\[\\[Alpha\]\] over the roots \\[Alpha\] of p\[\\[Alpha\]\] == 0.  Held symbolic form, used by the rational integrator's NaiveLogPart fallback when the logarithmic part cannot be expressed in closed-form real elementary functions. Differentiation threads through the body Function: D\[RootSum\[f1, Function\[t, body\]\], x\] == RootSum\[f1, Function\[t, D\[body, x\]\]\].

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= RootSum[Function[t, t^2 + t + 1], Function[t, t^3]]
Out[1]= RootSum[Function[t, t^2 + t + 1], Function[t, t^3]]

In[2]:= RootSum[#^2 + 1 &, # &]
Out[2]= RootSum[#1^2 + 1 &, #1 &]
```

```mathematica
In[1]:= Integrate[1/(x^3 - x + 1), x]
Out[1]= RootSum[1 + #1^3 - #1 &, Log[-#1 + x]/(-1 + 3 #1^2) &]
```

```mathematica
In[1]:= D[RootSum[Function[t, t^3 + t + 1], Function[t, Log[x - t]/t]], x]
Out[1]= RootSum[Function[t, t^3 + t + 1], Function[t, 1/(t (-t + x))]]
```

## Algorithm

root.c

```text
Mathematica-style symbolic Root and RootSum.  See root.h for the
representation contract.  Both heads are held forms — the C
```

builtins below intentionally return NULL for every input so the

```text
evaluator leaves the call exactly as the caller wrote it.  The
```

useful work happens elsewhere:

```text
  - src/deriv.c   — D[RootSum[f1, f2], x] threads through the body.
  - src/intrat.c  — NaiveLogPart constructs RootSum nodes when the
                    LRT log part has no closed-form real expression.
```

Splitting this module out of intrat.c keeps the held-symbolic machinery available to any future caller (e.g. Solve, Reduce, factorisation over algebraic extensions) that needs to name the roots of a polynomial without committing to a radical form.

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

## References

- Manuel Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed. (Springer, 2005).
- Source: [`src/root.c`](https://github.com/stblake/mathilda/blob/main/src/root.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)
- Tests: [`tests/test_intrat.c`](https://github.com/stblake/mathilda/blob/main/tests/test_intrat.c)

## Notes & additional examples

### Notes

`RootSum[f, form]` denotes the formal sum of `form[a]` over the roots `a` of
`f[a] == 0` and is kept as a held symbolic object. It is produced by the
rational integrator's `NaiveLogPart` fallback when the logarithmic part of an
integral cannot be written with real elementary functions; differentiation
threads through the `form` function.
