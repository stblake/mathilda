# ReplaceAt

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ReplaceAt[expr, rules, n]
    transforms expr by replacing the n-th element using rules.
ReplaceAt[expr, rules, {i, j, ...}]
    replaces the part of expr at position {i, j, ...}.
ReplaceAt[expr, rules, {{i1, j1, ...}, {i2, j2, ...}, ...}]
    replaces parts at several positions.

Rules may be a single Rule/RuleDelayed or a list of them; rules are tried
in order and the first match wins. Negative indices count from the end;
0 targets the head. All and Span specifications are supported. Repeated
positions cause rules to be applied repeatedly to that part.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ReplaceAt[{a, a, a, a}, a -> xx, 2]
Out[1]= {a, xx, a, a}

In[2]:= ReplaceAt[{a, a, a, a}, a -> xx, {{1}, {4}}]
Out[2]= {xx, a, a, xx}

In[3]:= ReplaceAt[{{a, a}, {a, a}}, a -> xx, {2, 1}]
Out[3]= {{a, a}, {xx, a}}

In[4]:= ReplaceAt[{a, a, a, a}, a -> xx, -2]
Out[4]= {a, a, xx, a}

In[5]:= ReplaceAt[{{a, a, a}, {a, a, a}}, a -> xx, {-1, -2}]
Out[5]= {{a, a, a}, {a, xx, a}}

In[6]:= ReplaceAt[{1, 2, 3, 4}, x_ :> 2 x - 1, {{2}, {4}}]
Out[6]= {1, 3, 3, 7}

In[7]:= ReplaceAt[{a, b, c, d}, {a -> xx, _ -> yy}, {{1}, {2}, {4}}]
Out[7]= {xx, yy, c, yy}

In[8]:= ReplaceAt[{{a, a}, {a, a}}, a -> xx, {All, 2}]
Out[8]= {{a, xx}, {a, xx}}
```

## Implementation notes

**Algorithm.** `builtin_replace_at` (`src/replace.c`) applies rules at one or more explicit *positions* rather than by structural matching everywhere. It parses the rule(s) with `parse_replace_rules` into a `ReplaceRule[]`, then disambiguates the position argument: a non-empty `List` whose first element is itself a `List` is a list of paths (applied sequentially, repeated positions re-apply the rules); otherwise it is a single path (a bare index or a `List` of indices). Navigation is `replaceat_at_path`, which consumes one index per level: index `0` descends into the head, a positive/negative integer selects an argument (negative counts from the end), `All` recurses into every argument, and a `Span[start, stop, step]` walks a strided slice. When the path is exhausted at a node, the rules are matched against *that node only* via the same `match`/`replace_bindings` machinery used by `Replace`/`ReplaceAll`; the first matching rule's bound replacement is substituted. Sub-trees off the targeted path are deep-copied unchanged.

**Data structures.** `ReplaceRule[]` (borrowed pattern/replacement pointers); the path is an `Expr**` slice advanced by pointer arithmetic (`path + 1`, `plen - 1`) as the recursion descends.

- `Protected`.
- `rules` may be a single `Rule` (`->`), `RuleDelayed` (`:>`), or a list of such rules. The rules are tried in order; the first one that applies wins. If no rule matches at a targeted position, the part is left unchanged.
- For `RuleDelayed`, the right-hand side is evaluated separately for each match after substituting bound pattern variables.
- Negative integer indices count from the end. The literal index `0` targets the head of an expression.
- Path components may be integers, the symbol `All` (selects every child at that level), or `Span` expressions such as `i ;; j` or `i ;; j ;; k`.
- Works on expressions with any head (not just `List`); after substitution the evaluator re-applies canonical ordering for `Orderless` heads such as `Plus` and `Times`.
- The position list uses the same form as is returned by `Position`. `ReplaceAt[expr, rules, {}]` applies the rules to the whole expression.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/replace.c`](https://github.com/stblake/mathilda/blob/main/src/replace.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
