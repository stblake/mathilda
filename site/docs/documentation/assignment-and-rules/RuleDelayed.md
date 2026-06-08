# RuleDelayed

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
lhs :> rhs or RuleDelayed[lhs, rhs]
    represents a delayed rewrite rule: rhs is held and evaluated only
    each time the rule fires, after the pattern bindings on lhs are
    substituted into rhs.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`RuleDelayed[lhs, rhs]` (`:>`) is a passive delayed-rewrite object with no builtin handler. Its `ATTR_HOLDREST | ATTR_SEQUENCEHOLD | ATTR_PROTECTED` attributes hold `rhs` unevaluated at construction; the rule engine (`is_rule` in `src/replace.c`) detects the `RuleDelayed` head and, each time the rule fires, substitutes the fresh pattern bindings into the held `rhs` and only then evaluates it (see the `delayed` flag threaded through `ReplaceRule`/`ReplaceListState`). This is the difference from `Rule`, whose RHS is evaluated once up front.

**Attributes:** `HoldRest`, `Protected`, `SequenceHold`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/replace.c`](https://github.com/stblake/mathilda/blob/main/src/replace.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= {1, 2, 3} /. n_ :> n^2
Out[1]= {1, 4, 9}

In[2]:= FullForm[a :> b]
Out[2]= RuleDelayed[a, b]
```

### Notes

`a :> b` is shorthand for `RuleDelayed[a, b]`. The right-hand side is held and evaluated separately for each match, after the pattern bindings are substituted in.
