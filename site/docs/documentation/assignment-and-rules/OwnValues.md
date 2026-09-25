# OwnValues

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`OwnValues[s] gives a list of own-value rules for s.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= f[x_] := x^2; DownValues[f]
Out[1]= {HoldPattern[f[x_]] :> x^2}

In[2]:= a = 5; OwnValues[a]
Out[2]= {HoldPattern[a] :> 5}
```

### Applications (3)

```mathematica
In[3]:= OwnValues[a]
Out[3]= {}

In[4]:= a = 7
Out[4]= 7

In[5]:= OwnValues[a]
Out[5]= {7 -> 7}
```

## Implementation notes

`builtin_own_values` (1-arg, symbol only) calls `symtab_get_own_values(sym)` to retrieve the symbol's `Rule*` list — the immediate `x = value` assignments — and hands it to the shared `rules_to_list` helper, which deep-copies each `pattern`/`replacement` pair into a `Rule[lhs, rhs]` node and gathers them into a `List`. The list is returned unevaluated so the stored values are reported verbatim rather than re-evaluated.

- Each rule is returned as `HoldPattern[lhs] :> rhs` (a `RuleDelayed` whose
  left-hand side is wrapped in `HoldPattern`), matching Mathematica. The
  `HoldPattern` keeps the stored pattern from re-matching its own definition, and
  the `RuleDelayed` keeps the right-hand side unevaluated, so the returned list is
  **inert**: querying `DownValues[f]` for a recursive `f` does not re-trigger the
  recursion, and a right-hand side such as `RandomInteger[...]` is shown held
  rather than evaluated afresh on each query.
- Returns `{}` for a symbol that carries no rules of the requested kind.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [DownValues](../../assignment-and-rules/DownValues/), [RuleDelayed](../../assignment-and-rules/RuleDelayed/), [HoldPattern](../../pattern-matching/HoldPattern/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_unset.c`](https://github.com/stblake/mathilda/blob/main/tests/test_unset.c)

## Notes & additional examples

### Notes

`OwnValues[s]` returns the direct value rules created by `s = ...`. An undefined symbol has an empty list; assigning `a = 7` stores a single own-value rule.
