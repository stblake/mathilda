# Blank

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
_ or Blank[] represents any single expression.
_h or Blank[h] represents any single expression with head h.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`Blank` is a pattern object head, not an evaluator builtin: `_` parses to `Blank[]` and `_h` to `Blank[h]`, and `x_` to `Pattern[x, Blank[...]]`. The matcher recognises it in `src/match.c` via `is_blank`, which reports whether a pattern node is `Blank[]` (matches any single expression — `head_out = NULL`) or `Blank[h]` (constrains the head). When a head is present, `match_internal` checks it with `get_expr_head_borrowed` + `expr_eq`, with a special case mapping atomic expression types to their symbolic heads (`Integer`, `Real`, `Symbol`, `String`). A successful single match consumes exactly one argument and continues via `call_parent`. `Blank` and its sequence relatives are also among the pattern heads the matcher treats structurally (e.g. excluded from `OneIdentity`/literal-head handling).

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/match.c`](https://github.com/stblake/mathilda/blob/main/src/match.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= MatchQ[5, _]
Out[1]= True

In[2]:= MatchQ[g[a], _g]
Out[2]= True

In[3]:= FullForm[x_]
Out[3]= Pattern[x, Blank[]]
```

### Notes

`_` (`Blank[]`) matches any single expression; `_h` (`Blank[h]`) matches any single expression whose head is `h` (Out[2]). A named blank `x_` is `Pattern[x, Blank[]]` (Out[3]), binding the match to `x`.
