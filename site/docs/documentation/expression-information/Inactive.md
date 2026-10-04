# Inactive

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Inactive[f] represents f with evaluation of its own rules suppressed, so Inactive[f][args] stays unevaluated (its arguments still evaluate).  Used to hold an integral inert, e.g. Inactive[Integrate][g, x]; Activate reverses it.`**

## Examples (6)

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

### Applications (4)

The head is held inert, so the sum is not taken

```mathematica
In[3]:= Inactive[Plus][2, 3]
Out[3]= Inactive[Plus][2, 3]
```

But the arguments still evaluate

```mathematica
In[4]:= Inactive[Plus][1 + 1, 3]
Out[4]= Inactive[Plus][2, 3]
```

A symbolic, inert integral

```mathematica
In[5]:= Inactive[Integrate][x, x]
Out[5]= Inactive[Integrate][x, x]
```

Activate turns it back on

```mathematica
In[6]:= Activate[Inactive[Integrate][x, x]]
Out[6]= 1/2 x^2
```

## Implementation notes

**What it is.** `Inactive[h]` is an inert wrapper for a head. It has no C evaluation
handler; it is registered in `core_init` (`src/core.c`) only so the symbol is official and
`ATTR_PROTECTED`. Inertness is automatic: `Inactive[h][args]` has the *compound* head
`Inactive[h]`, which carries no DownValues and no builtin, so the evaluator evaluates the
arguments normally but never fires `h`'s rules — the expression stays a symbolic fixed
point (the same mechanism as `Derivative[n][f][x]`).

**Interaction.** Only the head is held inert; the arguments still evaluate, so
`Inactive[Plus][1 + 1, 3]` first reduces to `Inactive[Plus][2, 3]`. `D` recognises inert
heads — it applies the fundamental theorem of calculus to `Inactive[Integrate]` (see
`src/calculus/deriv.c`) — and `Activate` reverses the wrapper by rewriting `Inactive[h] ->
h` throughout and re-evaluating. The interned pointer `SYM_Inactive` is the handle the
DSolve and derivative machinery use to build and detect inert integrals.

**Complexity / limits.** No cost of its own — it is a structural marker, matched by a
pointer comparison on the head. Its single subtlety is the head/argument asymmetry above:
the head's rules are suppressed, the arguments' are not.

- `Protected`.  Inertness is automatic: `Inactive[f][args]` is an unevaluated fixed point of a compound head (the head `Inactive[f]` carries no rule), analogous to `Derivative[n][f][x]`.
- Chief use: holding an integral inert.  `Inactive[Integrate][g, x]` never runs the integration cascade, so a non-elementary integrand (which would otherwise be expensive or non-terminating) is held instantly.
- **`D` applies the fundamental theorem of calculus** to an inactive integral WITHOUT evaluating it: `D[Inactive[Integrate][f, u], u]` = `f`.  (A different differentiation variable falls to the ordinary rules; an integrand free of it gives `0`.)  This lets an inert first integral verify by the implicit-function rule with no integration cost.
- `Activate` reverses it.

**Attributes:** `Protected`.

## References

**See also:** [D](../../calculus/D/), [Activate](../../expression-information/Activate/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_deriv.c`](https://github.com/stblake/mathilda/blob/main/tests/test_deriv.c)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)
- Tests: [`tests/test_dsolve_m58_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m58_stress.c)
- Tests: [`tests/test_dsolve_m61_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m61_stress.c)

## Notes & additional examples

### Notes

`Inactive[h]` wraps a head so that its own rules are suppressed: `Inactive[h][args]` has
the compound head `Inactive[h]`, which has no rules of its own, so the call stays a
symbolic fixed point. It is how a computation like an integral is kept inert for display
or later manipulation.

Only the head is held — the arguments still evaluate, so `Inactive[Plus][1 + 1, 3]` first
reduces to `Inactive[Plus][2, 3]`. `Activate` reverses the wrapper, and `D` already knows
how to differentiate an inert `Inactive[Integrate]`.
