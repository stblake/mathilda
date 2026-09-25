# Inactive

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Inactive[f] represents f with evaluation of its own rules suppressed, so Inactive[f][args] stays unevaluated (its arguments still evaluate).  Used to hold an integral inert, e.g. Inactive[Integrate][g, x]; Activate reverses it.`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

Held; no integration

```mathematica
In[1]:= Inactive[Integrate][1/Sqrt[y Log[y] + 3], y]
Out[1]= Inactive[Integrate][1/Sqrt[3 + y Log[y]], y]
```

```mathematica
In[2]:= D[Inactive[Integrate][1/p[y], y], y]
Out[2]= 1/p[y]
```

## Implementation notes

- `Protected`.  Inertness is automatic: `Inactive[f][args]` is an unevaluated fixed point of a compound head (the head `Inactive[f]` carries no rule), analogous to `Derivative[n][f][x]`.
- Chief use: holding an integral inert.  `Inactive[Integrate][g, x]` never runs the integration cascade, so a non-elementary integrand (which would otherwise be expensive or non-terminating) is held instantly.
- **`D` applies the fundamental theorem of calculus** to an inactive integral WITHOUT evaluating it: `D[Inactive[Integrate][f, u], u]` = `f`.  (A different differentiation variable falls to the ordinary rules; an integrand free of it gives `0`.)  This lets an inert first integral verify by the implicit-function rule with no integration cost.
- `Activate` reverses it.

**Attributes:** `Protected`.

## References

**See also:** [D](../../calculus/D/), [Activate](../../expression-information/Activate/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)
- Tests: [`tests/test_dsolve_m58_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m58_stress.c)
- Tests: [`tests/test_inactive.c`](https://github.com/stblake/mathilda/blob/main/tests/test_inactive.c)
