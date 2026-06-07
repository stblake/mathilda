# ReplaceAll

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
expr /. rules or ReplaceAll[expr, rules]
    traverses expr top-down and applies the first matching rule at each
    subexpression. A matched subexpression is replaced and NOT recursed
    into further -- ReplaceAll is a single pass, not a fixed point.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= {x, x^2, y, z} /. x -> 1
Out[1]= {1, 1, y, z}

In[2]:= Sin[x] /. Sin -> Cos
Out[2]= Cos[x]

In[3]:= {1, 3, 2, x, 6, Pi} /. _?PrimeQ -> "prime"
Out[3]= {1, "prime", "prime", x, 6, Pi}

In[4]:= {f[2], f[x, y], h[], f[]} /. f[x__] -> "OK"
Out[4]= {"OK", "OK", h[], f[]}
```

## Implementation notes

**Algorithm.** `builtin_replace_all` (`/.`, in `src/replace.c`) drives `apply_replace_all_nested`, which first collects the rule(s) into a flat `ReplaceRule[]` (pattern, replacement, delayed flag). A `List` of rules whose first element is *not* itself a `Rule`/`RuleDelayed` is treated as a list of alternative rule-sets and threaded — one result per sub-list. The core traversal is `do_replace_all`: at each node it tries every rule in order via `match(e, rule.pattern, env)` (the structural matcher in `src/match.c`); on the first match it builds the result with `replace_bindings(rule.replacement, env)` and **returns immediately without descending into the matched subexpression** — this is the defining top-down, non-re-entrant `ReplaceAll` semantic. If no rule matches the node, it recurses into the head and every argument and rebuilds the function node. Rule-creation timing (immediate `Rule` vs delayed `RuleDelayed`) is honoured by `replace_bindings`/the surrounding evaluator rather than re-applied here.

**Data structures.** `ReplaceRule` array holding borrowed pointers into the user's rule expressions; a fresh `MatchEnv` (from `env_new`) per match attempt holds the pattern-variable bindings consumed by `replace_bindings`.

- `Protected`.
- Evaluates the entire expression top-down. The first rule that applies to a particular part is used; no further rules are tried on that part or on any of its subparts.
- Applies a rule only once to an expression.
- Returns `expr` unmodified if no rules apply.
- Maps across lists of rules appropriately.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/replace.c`](https://github.com/stblake/mathilda/blob/main/src/replace.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= x + y /. x -> 2
Out[1]= 2 + y
```

```mathematica
In[1]:= {x, y, z} /. {x -> 1, z -> 3}
Out[1]= {1, y, 3}
```

```mathematica
In[1]:= x^2 + x /. x -> a + 1
Out[1]= 1 + a + (1 + a)^2
```

```mathematica
In[1]:= f[1, 2] /. f[a_, b_] -> a + b
Out[1]= 3
```

### Notes

`expr /. rules` is the shorthand for `ReplaceAll[expr, rules]`. It traverses
`expr` top-down and rewrites each subexpression that matches a rule's left-hand
side. A single rule or a list of rules may be supplied; with a list, the first
matching rule wins at each position. Replacement is not re-applied to the result
of a match (that is `ReplaceRepeated`, `//.`). Patterns with named blanks
(`f[a_, b_] -> a + b`) bind variables for use on the right-hand side.
