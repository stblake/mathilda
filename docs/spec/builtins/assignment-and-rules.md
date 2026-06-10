# Assignment and Rules

## Set (=), SetDelayed (:=)
- `lhs = rhs`: `Set`. Immediate evaluation of `rhs`.
- `lhs := rhs`: `SetDelayed`. Delayed evaluation of `rhs`.
- `lhs := rhs /; condition`: When the RHS of `SetDelayed` is a `Condition`, it is automatically moved to the LHS pattern. This makes `f[x_] := body /; test` equivalent to `f[x_] /; test := body`.

## Unset (=.)

- `Unset[lhs]` or `lhs =.`: Removes any rule whose left-hand side is `lhs`. A
  bare symbol clears its `OwnValue`; a function form (`f[args]`) clears the one
  matching `DownValue` on the head symbol, leaving sibling rules untouched.
- Rules are removed only when their left-hand sides are identical to `lhs` *up
  to renaming of bound pattern variables* — so `f[x_] =.` removes a rule defined
  as `f[y_] := ...`, while `fact[1] =.` leaves the separate `fact[n_]` rule in
  place.

**Features**:
- `=.` is a low-precedence postfix operator (precedence 40, like `Set`), so it
  captures the whole preceding expression: `a b =.` parses as `Unset[a b]`. The
  guard against a trailing digit keeps `k =.5` parsing as `Set[k, 0.5]`.
- `Unset` has attributes `{HoldFirst, Protected}`; it holds `lhs`, so the symbol
  (not its value) is operated on. `Protected`/`Locked` symbols are not affected.
- Always returns `Null`, whether or not a matching rule was found.

```mathematica
In[1]:= x = 5; x =.; x
Out[1]= x

In[2]:= f[x_] := x^2; f[x_] =.; f[3]
Out[2]= f[3]

In[3]:= fact[1] = 1; fact[n_] := n fact[n - 1]; fact[1] =.; fact[1]
Out[3]= fact[1]
```

## Replace
Applies a rule or list of rules to transform an expression.
- `Replace[expr, rules]`: Applies rules to the entire expression.
- `Replace[expr, rules, levelspec]`: Applies rules to parts specified by `levelspec`.
- `Replace[expr, rules, levelspec, Heads -> True]`: Includes heads of expressions and their parts.

**Features**:
- `Protected`.
- Rules must be of the form `lhs -> rhs` (`Rule`) or `lhs :> rhs` (`RuleDelayed`).
- Tries rules in order. The first one that matches is applied.
- If rules are given in nested lists, `Replace` is mapped onto the inner lists.
- Standard level specifications (`n`, `Infinity`, `All`, `{n}`, `{n1, n2}`) are fully supported with the default being `{0}` (the whole expression).
- Expressions at deeper levels in a subexpression are matched first.
- Replaces parts even inside `Hold` or related wrappers.

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

## ReplacePart
Replaces parts of an expression at specified positions.
- `ReplacePart[expr, i -> new]`: Replaces the part at position `i`.
- `ReplacePart[expr, {i, j, ...} -> new]`: Replaces the part at a nested position.
- `ReplacePart[expr, {{i1, ...} -> new1, {i2, ...} -> new2, ...}]`: Multiple replacements.
- `ReplacePart[expr, pattern -> new]`: Replaces parts at all positions matching `pattern` (e.g., `_`, `Except[...]`).

**Features**:
- `Protected`.
- Supports negative indices: `-1` is the last element, `-2` is second-to-last, etc.
- Supports mixed positive/negative indices in nested position specs (e.g., `{2, -1}`).
- Default `Heads -> False`: pattern-based position specs do not match the head (index 0). Use an explicit `0 -> new` rule to replace the head.
- Supports `Rule` (`->`) and `RuleDelayed` (`:>`).

```mathematica
In[1]:= ReplacePart[{a, b, c, d}, 2 -> x]
Out[1]= {a, x, c, d}

In[2]:= ReplacePart[{a, b, c, d, e, f, g}, -3 -> xxx]
Out[2]= {a, b, c, d, xxx, f, g}

In[3]:= ReplacePart[{a, b, c, d, e, f, g}, Except[1 | 3 | 5] -> xxx]
Out[3]= {a, xxx, c, xxx, e, xxx, xxx}

In[4]:= ReplacePart[f[x], 0 -> g]
Out[4]= g[x]
```

## ReplaceAt
Replaces parts of an expression at specified positions using replacement rules. Unlike `ReplaceAll` (which traverses the whole tree), `ReplaceAt` only attempts the rules at the part(s) named by the position specification.

- `ReplaceAt[expr, rules, n]`: replaces the `n`-th element using `rules`.
- `ReplaceAt[expr, rules, {i, j, ...}]`: replaces the part of `expr` at position `{i, j, ...}` (equivalently `expr[[i, j, ...]]`).
- `ReplaceAt[expr, rules, {{i1, j1, ...}, {i2, j2, ...}, ...}]`: replaces the parts at each of the listed positions. If the same position appears more than once, the rules are applied repeatedly to that part.

**Features**:
- `Protected`.
- `rules` may be a single `Rule` (`->`), `RuleDelayed` (`:>`), or a list of such rules. The rules are tried in order; the first one that applies wins. If no rule matches at a targeted position, the part is left unchanged.
- For `RuleDelayed`, the right-hand side is evaluated separately for each match after substituting bound pattern variables.
- Negative integer indices count from the end. The literal index `0` targets the head of an expression.
- Path components may be integers, the symbol `All` (selects every child at that level), or `Span` expressions such as `i ;; j` or `i ;; j ;; k`.
- Works on expressions with any head (not just `List`); after substitution the evaluator re-applies canonical ordering for `Orderless` heads such as `Plus` and `Times`.
- The position list uses the same form as is returned by `Position`. `ReplaceAt[expr, rules, {}]` applies the rules to the whole expression.

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

## ReplaceAll (/.)
Applies a rule or list of rules to transform each subpart of an expression.
- `expr /. lhs -> rhs` or `expr /. {rules}`

**Features**:
- `Protected`.
- Evaluates the entire expression top-down. The first rule that applies to a particular part is used; no further rules are tried on that part or on any of its subparts.
- Applies a rule only once to an expression.
- Returns `expr` unmodified if no rules apply.
- Maps across lists of rules appropriately.

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

## CompoundExpression (;)
- `expr1; expr2; ...`: Evaluates a sequence of expressions, returning the last one.

## ClearAll, Remove

- `ClearAll[s1, s2, ...]`: Clears all values (OwnValues and DownValues),
  attributes and usage messages for the named symbols, leaving the symbols
  themselves in the symbol table. `ClearAll[{s1, s2, ...}]` accepts a list of
  specs.
- `Remove[s1, s2, ...]`: Removes the named symbols completely, deleting their
  definitions from the symbol table. A later reference recreates a fresh, empty
  symbol of the same name. `Remove[{s1, s2, ...}]` accepts a list of specs.

**Features**:
- `ClearAll` has attributes `{HoldAll, Protected}`; `Remove` has
  `{HoldAll, Locked, Protected}`. Both hold their arguments, so they operate on
  the symbol, not its current value.
- Neither affects symbols with the attribute `Locked` or `Protected`. This is
  what prevents `Remove`/`ClearAll` from ever deleting or wiping a built-in.
- `ClearAll`, unlike `Clear`, also removes attributes and the usage message.
- Both return `Null`.

```mathematica
In[1]:= f[x_] := x^2; SetAttributes[f, Listable]; Attributes[f]
Out[1]= {Listable}

In[2]:= ClearAll[f]; {Attributes[f], DownValues[f]}
Out[2]= {{}, {}}

In[3]:= x = 2; Remove[x]; x
Out[3]= x
```

## Protect, Unprotect

- `Protect[s1, s2, ...]`: Sets the attribute `Protected` for the named symbols.
- `Unprotect[s1, s2, ...]`: Removes the attribute `Protected` from the named
  symbols.
- Both accept a list of specs (`Protect[{s1, s2, ...}]`) and return the list of
  names (as strings) whose protection state actually changed — `{}` when nothing
  changed.

**Features**:
- Both have attributes `{HoldAll, Protected}` and hold their arguments.
- Neither affects symbols with the attribute `Locked`.
- The typical sequence for adding rules to an existing symbol is
  `Unprotect[f]; definition; Protect[f]`.

```mathematica
In[1]:= f[x_] := x^2; Protect[f]
Out[1]= {"f"}

In[2]:= f[x_, y_] := x + y
SetDelayed::wrsym: Symbol f is Protected.

In[3]:= Unprotect[f]
Out[3]= {"f"}
```

