# ReplaceAt

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ReplaceAt[expr, rules, n]`**

transforms expr by replacing the n-th element using rules.

**`ReplaceAt[expr, rules, {i, j, ...}]`**

replaces the part of expr at position {i, j, ...}.

**`ReplaceAt[expr, rules, {{i1, j1, ...}, {i2, j2, ...}, ...}]`**

replaces parts at several positions.

<details>
<summary>Notes</summary>

Rules may be a single Rule/RuleDelayed or a list of them; rules are tried in order and the first match wins. Negative indices count from the end; 0 targets the head. All and Span specifications are supported. On an association a position is a key, Key\[k\], or a positional index over the entries, and the rules are tried against the value. Repeated positions cause rules to be applied repeatedly to that part. ReplaceAt\[expr, rules, {}\] is an empty list of positions and replaces nothing, while {{}} is the position of expr itself. A position that does not exist leaves ReplaceAt unevaluated.

</details>

## Examples (22)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (17)

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

In[9]:= ReplaceAt[{{a, b}, {c, d}, e}, x_ :> f[x], 2]
Out[9]= {{a, b}, f[{c, d}], e}

In[10]:= ReplaceAt[{{a, b}, {c, d}, e}, x_ :> f[x], -1]
Out[10]= {{a, b}, {c, d}, f[e]}

In[11]:= ReplaceAt[{{a, b}, {c, d}, e}, x_ :> f[x], {2, 1}]
Out[11]= {{a, b}, {f[c], d}, e}

In[12]:= ReplaceAt[{{a, b}, {c, d}, e}, x_ :> f[x], {{1}, {3}}]
Out[12]= {f[{a, b}], {c, d}, f[e]}

In[13]:= ReplaceAt[{{a, b}, {c, d}, e}, x_ :> f[x], {{1, 2}, {2, 2}, {3}}]
Out[13]= {{a, f[b]}, {c, f[d]}, f[e]}

In[14]:= ReplaceAt[{a, a, a, a, a}, a -> xx, 2 ;; 4]
Out[14]= {a, xx, xx, xx, a}

In[15]:= ReplaceAt[a + b + c + d, _ -> x, 2]
Out[15]= a + c + d + x

In[16]:= ReplaceAt[x^2 + y^2, _ -> z, {{1, 1}, {2, 1}}]
Out[16]= 2 z^2

In[17]:= ReplaceAt[{a, b, c}, _ -> f, 0]
Out[17]= f[a, b, c]
```

### Applications (5)

```mathematica
In[18]:= ReplaceAt[{a, b, c, d}, x_ -> X, 2]
Out[18]= {a, X, c, d}

In[19]:= ReplaceAt[{a, b, c, d}, x_ -> X, -1]
Out[19]= {a, b, c, X}

In[20]:= ReplaceAt[{{a, b}, {c, d}}, x_ -> X, {2, 1}]
Out[20]= {{a, b}, {X, d}}

In[21]:= ReplaceAt[1 + x + x^2 + x^3, e_ :> D[e, x], {2}]
Out[21]= 2 + x^2 + x^3

In[22]:= ReplaceAt[{1, 2, 3, 4, 5}, n_ :> n^2, {2 ;; 4}]
Out[22]= {1, 4, 9, 16, 5}
```

## Implementation notes

**Algorithm.** `builtin_replace_at` (`src/replace.c`) applies rules at one or more explicit *positions* rather than by structural matching everywhere. It parses the rule(s) with `parse_replace_rules` into a `ReplaceRule[]`, then disambiguates the position argument: a non-empty `List` whose first element is itself a `List` is a list of paths (applied sequentially, repeated positions re-apply the rules); otherwise it is a single path (a bare index or a `List` of indices). Navigation is `replaceat_at_path`, which consumes one index per level: index `0` descends into the head, a positive/negative integer selects an argument (negative counts from the end), `All` recurses into every argument, and a `Span[start, stop, step]` walks a strided slice. When the path is exhausted at a node, the rules are matched against *that node only* via the same `match`/`replace_bindings` machinery used by `Replace`/`ReplaceAll`; the first matching rule's bound replacement is substituted. Sub-trees off the targeted path are deep-copied unchanged.

**Data structures.** `ReplaceRule[]` (borrowed pattern/replacement pointers); the path is an `Expr**` slice advanced by pointer arithmetic (`path + 1`, `plen - 1`) as the recursion descends.

- `Protected`.
- `rules` may be a single `Rule` (`->`), `RuleDelayed` (`:>`), or a list of such rules. The rules are tried in order; the first one that applies wins. If no rule matches at a targeted position, the part is left unchanged.
- For `RuleDelayed`, the right-hand side is evaluated separately for each match after substituting bound pattern variables.
- Negative integer indices count from the end. The literal index `0` targets the head of an expression.
- Path components may be integers, the symbol `All` (selects every child at that level), or `Span` expressions such as `i ;; j` or `i ;; j ;; k` (including `UpTo[n]` bounds).
- Works on expressions with any head (not just `List`); after substitution the evaluator re-applies canonical ordering for `Orderless` heads such as `Plus` and `Times`.
- On an `Association` a position is a key, `Key[k]`, or a positional index over the entries, and the rules are tried against the *value*; `All`, `Span` and `0` work there too.
- The position list uses the same form as is returned by `Position`. `ReplaceAt[expr, rules, {}]` is an **empty list of positions** and replaces nothing; the position of the whole expression is the empty path, `{{}}`.
- A position that does not exist — an out-of-range index, an absent key, a malformed `Span`, or a path that runs into an atom — leaves `ReplaceAt` unevaluated, as `Part` does for `{a, b, c}[[99]]`.
- Position resolution is shared with [`MapAt`](../data-structures/MapAt.md) via one walker (`expr_apply_at_path`, `src/part.c`), so the two agree on every position spec by construction.

**Attributes:** `Protected`.

## References

**See also:** [ReplaceAll](../../assignment-and-rules/ReplaceAll/), [Rule](../../assignment-and-rules/Rule/), [RuleDelayed](../../assignment-and-rules/RuleDelayed/), [Span](../../structural-manipulation/Span/), [List](../../other-advanced/List/), [Orderless](../../expression-information/Orderless/), [Plus](../../arithmetic/Plus/), [Times](../../arithmetic/Times/)

- Source: [`src/replace.c`](https://github.com/stblake/mathilda/blob/main/src/replace.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_replaceat.c`](https://github.com/stblake/mathilda/blob/main/tests/test_replaceat.c)

## Notes & additional examples

### Notes

`ReplaceAt` applies rules only at explicitly named positions, leaving the rest
of the expression untouched. Positions may be a single index, a part list
`{i, j, ...}`, a list of positions, `All`, or a `Span`. Negative indices count
from the end and `0` targets the head. Because the rule sees the *part* as a
whole expression, position-targeted rewrites can do real work: applying
`e_ :> D[e, x]` at part `{2}` of `1 + x + x^2 + x^3` differentiates just that
one summand (here the `x` term, whose derivative `1` merges into the leading
constant), while `n_ :> n^2` over the span `2 ;; 4` squares a contiguous slice.
