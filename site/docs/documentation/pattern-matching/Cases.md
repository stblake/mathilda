# Cases

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Cases[{e1, e2, ...}, pattern] gives a list of the ei that match the pattern.
Cases[{e1, ...}, pattern -> rhs] gives a list of the values of rhs corresponding to the ei that match the pattern.
Cases[expr, pattern, levelspec] gives a list of all parts of expr on levels specified by levelspec that match the pattern.
Cases[expr, pattern -> rhs, levelspec] gives the values of rhs that match the pattern.
Cases[expr, pattern, levelspec, n] gives the first n parts in expr that match the pattern.
Cases[pattern] represents an operator form of Cases that can be applied to an expression.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_cases` (`src/patterns.c`) collects subexpressions matching a pattern. It parses the level-spec (default `{1,1}` — immediate elements only; integer `n` → `1..n`, negative `n` → exactly level `n` from the bottom, `{m,n}`, `All`/`Infinity`), a `Heads -> True|False` option, and an optional `n` result limit. If the second argument is `pat -> rhs` / `pat :> rhs` it splits into a pattern plus a replacement applied to each hit. The traversal `do_cases_at_level` is depth-first **pre-order** (it recurses into head — when `heads` — and arguments first, then tests the node): at every in-range node it calls `match(e, pattern, env)` from `src/match.c`; on success it appends either a deep copy of the node, or, for a rule, `replace_bindings(rhs, env)` (evaluated when the rule is `RuleDelayed`). Level membership for negative specs is decided against `get_expr_depth_patterns`. Collection stops once `max_results` hits are gathered. `Cases[pat]` with one argument returns an operator form `Function[Cases[#1, pat]]`.

**Data structures.** A growable `Expr**` results buffer (doubled on demand) wrapped into a `List` at the end; one `MatchEnv` per node tested.

- Uses standard level specifications.
- Defaults to level `{1}`.
- Option `Heads -> True` includes heads of expressions.
- Operates in depth-first post-order traversal.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/patterns.c`](https://github.com/stblake/mathilda/blob/main/src/patterns.c)
- Specification: [`docs/spec/builtins/pattern-matching.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/pattern-matching.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Cases[{1, a, 2, b, 3}, _Integer]
Out[1]= {1, 2, 3}
```

```mathematica
In[1]:= Cases[{1, 2, 3, 4}, x_ /; x > 2]
Out[1]= {3, 4}
```

```mathematica
In[1]:= Cases[{f[1], g[2], f[3]}, f[x_] -> x]
Out[1]= {1, 3}
```

```mathematica
In[1]:= Cases[{{1, 2}, {3, 4}}, _Integer, 2]
Out[1]= {1, 2, 3, 4}
```

### Notes

`Cases[list, pattern]` returns the elements that match `pattern`. With a
`pattern -> rhs` rule it instead returns the transformed values (`f[x_] -> x`
extracts the argument of every `f`). Conditional patterns (`x_ /; x > 2`) filter
on a predicate. A trailing level specification (here `2`) descends into nested
lists and collects matches from those deeper levels — without it `Cases` only
inspects level 1.
