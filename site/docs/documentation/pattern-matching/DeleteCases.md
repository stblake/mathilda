# DeleteCases

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
DeleteCases[expr, pattern] removes all elements of expr that match pattern.
DeleteCases[expr, pattern, levelspec] removes all parts of expr on levels specified by levelspec that match pattern.
DeleteCases[expr, pattern, levelspec, n] removes the first n parts of expr that match pattern.
DeleteCases[pattern] represents an operator form of DeleteCases that can be applied to an expression.
The default levelspec is {1}. With Heads -> True, the heads of expressions are also tested; deleting a head is equivalent to applying FlattenAt at that location.
DeleteCases traverses expr in depth-first post-order (leaves before roots).
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_delete_cases` (`src/patterns.c`) returns a copy of the input with every subexpression matching the pattern (within the level-spec) removed. Level-spec, `Heads` option, and an optional deletion-count budget `n` are parsed as in `Cases` (default level `{1,1}`). The worker `do_delete_cases_at_level` is a depth-first **post-order** rebuild ("leaves before roots"): it recursively transforms the head (only if `heads`) and each argument, dropping any child flagged for deletion and splicing any `Sequence[...]` a head-deletion produced, then rebuilds the node. After rebuilding, it tests the **original** node against the pattern via `match`; a match sets the parent's `*delete_me` flag so the node is dropped from its parent's argument list. A matching head (under `Heads -> True`) turns the call into `Sequence[args…]` (FlattenAt-style), which the enclosing loop splices outward. The `n` budget (`count_remaining`: `-1` unlimited, else decremented per deletion) caps the number removed.

**Data structures.** Per-node freshly-allocated `Expr**` argument buffer (grown as needed); a `bool* delete_me` out-parameter threads the removal decision up to the parent; one `MatchEnv` per match test.

- Uses standard level specifications, defaulting to level `{1}`.
- Option `Heads -> True` tests heads as well; deleting a head is equivalent to applying `FlattenAt` at that point, splicing the remaining arguments into the parent.
- Traverses `expr` in depth-first post-order (leaves before roots) so that the `n` budget is spent on deeper matches before shallower ones.
- The match test is applied to the original subexpression (not the version with children already deleted), matching Mathematica semantics.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/patterns.c`](https://github.com/stblake/mathilda/blob/main/src/patterns.c)
- Specification: [`docs/spec/builtins/pattern-matching.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/pattern-matching.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= DeleteCases[{1,2,3,2},2]
Out[1]= {1, 3}
```

```mathematica
In[1]:= DeleteCases[{1,a,2,b,3},_Integer]
Out[1]= {a, b}
```

```mathematica
In[1]:= DeleteCases[{{1,2},{a,b},{3,c}}, {_Integer, _Integer}]
Out[1]= {{a, b}, {3, c}}
```

```mathematica
In[1]:= DeleteCases[{a + b, c d, e^2, f}, _Power, Infinity]
Out[1]= {a + b, c d, f}
```

### Notes

Removes every element matching the pattern; combine with `_Integer` or `_?OddQ` to filter by head or predicate. The pattern can be structural — `{_Integer, _Integer}` deletes only the sublists that are pairs of integers. With a level specification such as `Infinity`, `DeleteCases` descends into subexpressions and removes matching parts at any depth (here every `Power` subterm), traversing in depth-first post-order. The default level is `{1}` (top-level elements).
