# ReplaceList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ReplaceList[expr, rules] attempts to transform the entire expression expr by applying a rule or list of rules in all possible ways, and returns a list of the results obtained.
ReplaceList[expr, rules, n] gives a list of at most n results.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_replacelist` returns *all* ways a rule (or rule list) can match `expr`, rather than the first. It flattens the rule argument into a `ReplaceRule[]` array (recording the `delayed` flag for `RuleDelayed`), then for each rule runs the pattern matcher with a callback installed on the `MatchEnv`. Instead of stopping at the first successful binding, the matcher's backtracking enumerates every distinct binding (notably every partition of a `BlankSequence`/`BlankNullSequence`), and `replacelist_callback` materialises each one: it builds the replacement via `replace_bindings`, additionally `eval_and_free`-ing it for delayed rules, and appends it to a growable `ReplaceListState.results` buffer. An optional third argument caps the number of results (`state.limit`), checked both before invoking the matcher's callback and between rules. The accumulated results are wrapped in a `List`.

**Data structures.** `ReplaceRule { pattern, replacement, delayed }` and `ReplaceListState { results, count, cap, limit, replacement, delayed }` — the matcher signals each match through `env->callback`/`env->callback_data`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/replace.c`](https://github.com/stblake/mathilda/blob/main/src/replace.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= ReplaceList[{a, b, c, d}, {x___, y___} :> {{x}, {y}}]
Out[1]= {{{}, {a, b, c, d}}, {{a}, {b, c, d}}, {{a, b}, {c, d}}, {{a, b, c}, {d}}, {{a, b, c, d}, {}}}

In[2]:= ReplaceList[a + b + c, x_ + y_ :> {x, y}]
Out[2]= {{a, b + c}, {b, a + c}, {c, a + b}, {a + b, c}, {a + c, b}, {b + c, a}}

In[3]:= ReplaceList[{a, b, c, d}, {x___, y___} :> {{x}, {y}}, 2]
Out[3]= {{{}, {a, b, c, d}}, {{a}, {b, c, d}}}
```

```mathematica
In[1]:= ReplaceList[{1, 2, 3, 4}, {x___, y_, z_, w___} /; y + z == 5 :> {y, z}]
Out[1]= {{2, 3}}
```

```mathematica
In[1]:= ReplaceList[{1, 2, 3, 4, 5}, {a___, b_, c___} /; b == 3 :> {{a}, {c}}]
Out[1]= {{{1, 2}, {4, 5}}}
```

### Notes

Unlike `Replace`, which returns only the first match, `ReplaceList` enumerates *every* way the rule can match — useful for sequence patterns (`__`, `___`) that admit multiple partitions. A third argument caps the number of results returned.
