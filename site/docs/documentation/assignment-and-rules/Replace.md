# Replace

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Replace[expr, rules]
    tries to match expr at the top level against rules and returns the
    rewritten form; if no rule matches, returns expr unchanged.
Replace[expr, rules, levelspec]
    applies rules only at the parts of expr specified by levelspec.
Matching tries each rule in order and uses the first that succeeds;
rules may be a single rule or a list of rules.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Replace[x^2, x^2 -> a + b]
Out[1]= a + b

In[2]:= Replace[1 + x^2, x^2 -> a + b]
Out[2]= 1 + x^2

In[3]:= Replace[x, {{x -> a}, {x -> b}}]
Out[3]= {a, b}

In[4]:= Replace[1 + x^2, x^2 -> a + b, {1}]
Out[4]= 1 + a + b
```

## Implementation notes

**Algorithm.** `builtin_replace` (`src/replace.c`) differs from `ReplaceAll` in that it applies rules only at the levels named by an optional level-spec (third argument). The handler parses the level-spec into `[min_l, max_l]` (an integer `n` → levels `1..n`; `{n}` → exactly `n`; `{m,n}` → `m..n`; `All`/`Infinity` map to wide bounds; negative bounds select by depth-from-bottom), plus a `Heads -> True|False` option. It then calls `apply_replace_nested` → `do_replace_at_level`, a bottom-up traversal: it first recurses into the head (only if `heads`) and all arguments, rebuilds the node, and *then*, if the current node's level lies in `[min_l, max_l]` (level computed against `get_expr_depth_replace` for negative specs), tries each `ReplaceRule` in order via `match`/`replace_bindings` and returns the first replacement. Because matching happens after the children are already transformed and a match short-circuits further rule attempts at that node, `Replace` rewrites each in-range position once. With no level-spec it defaults to level 0 only (the whole expression).

**Data structures.** `ReplaceRule[]` of borrowed pattern/replacement pointers; `MatchEnv` per attempt. A list of non-rule sub-lists is threaded into parallel results, as in `ReplaceAll`.

- `Protected`.
- Rules must be of the form `lhs -> rhs` (`Rule`) or `lhs :> rhs` (`RuleDelayed`).
- Tries rules in order. The first one that matches is applied.
- If rules are given in nested lists, `Replace` is mapped onto the inner lists.
- Standard level specifications (`n`, `Infinity`, `All`, `{n}`, `{n1, n2}`) are fully supported with the default being `{0}` (the whole expression).
- Expressions at deeper levels in a subexpression are matched first.
- Replaces parts even inside `Hold` or related wrappers.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/replace.c`](https://github.com/stblake/mathilda/blob/main/src/replace.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Replace[x^2, x^2 -> done]
Out[1]= done
```

```mathematica
In[1]:= Replace[x, {x -> 1, _ -> 0}]
Out[1]= 1

In[2]:= Replace[w, {x -> 1, _ -> 0}]
Out[2]= 0
```

```mathematica
In[1]:= Replace[a + b + c + d, x_ + y_ -> {x, y}]
Out[1]= {a, b + c + d}
```

```mathematica
In[1]:= Replace[{{1, 2}, {3, 4}}, x_Integer :> x^2, {2}]
Out[1]= {{1, 4}, {9, 16}}
```

```mathematica
In[1]:= Replace[{1, {2, {3, {4}}}}, x_Integer :> x^2, {-1}]
Out[1]= {1, {4, {9, {16}}}}
```

### Notes

`Replace[expr, rules]` matches `expr` only at the **top level** (unlike
`ReplaceAll`, which descends into every subexpression). The first matching
rule in the list wins; a fall-through `_ -> 0` rule acts as a default. A
`levelspec` third argument restricts where matching happens: `{2}` targets
parts two levels deep, and the negative spec `{-1}` targets the leaves, so
`x_Integer :> x^2` squares every atom in a nested list. Because matching
respects flat/orderless heads, `x_ + y_` against a four-term sum binds `x` to
the first term and `y` to the rest.
