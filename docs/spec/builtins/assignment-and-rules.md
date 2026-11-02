# Assignment and Rules

## Set (=), SetDelayed (:=)
- `lhs = rhs`: `Set`. Immediate evaluation of `rhs`.
- `lhs := rhs`: `SetDelayed`. Delayed evaluation of `rhs`.
- `lhs := rhs /; condition`: When the RHS of `SetDelayed` is a `Condition`, it is automatically moved to the LHS pattern. This makes `f[x_] := body /; test` equivalent to `f[x_] /; test := body`.
- **The left-hand side is evaluated**: its head (`h = g;
  h[x_] := ...` defines `g`, and so does a scoped local holding a computed
  symbol, `ch = ToExpression["cache" <> name]; ch[a_] := ...`), and any pattern
  argument that mentions a symbol with a value (`PD = PropagatorDenominator;
  f[a__PD] := ...` is a rule for `__PropagatorDenominator`). In any other
  pattern argument, every pattern-free subexpression is evaluated, as it would
  be outside a pattern, while the pattern nodes stay as written: `h[y_ - 4] :=
  ...` is stored as `h[-4 + y_]` and fires on `h[x - 4]`, and `x_/2` holds
  `1/2`. Pattern names, guards (`/;`), tests (`?`), `HoldPattern` bodies and
  the held arguments of a Hold-attributed head are not evaluated (v0.586,
  v0.589).
- **`HoldPattern[f[...]] := rhs`** is a rule for `f`, stored with its
  `HoldPattern` (also under a guard, `HoldPattern[f[x_]] /; c := ...`, or with
  one inside); `HoldPattern[s] = v` is an OwnValue of `s`. Before v0.591 the
  rule silently landed on `HoldPattern`.
- **List destructuring.** `{a, b} = {1, 2}` assigns element by element and
  `{a, b} = c` threads `c` to each target. Lists of different lengths are
  `Set::shape` ("Lists {a, b} and {1} are not the same shape.") and nothing is
  assigned; the value of the `Set` is the evaluated right side, and of a
  `SetDelayed`, `$Failed` (v0.6144; it previously stayed unevaluated in
  silence). A mismatch nested inside a List is reported the same way. Here the
  reference still binds the siblings matched before the mismatch; Mathilda binds
  none of them.
- **Protected targets.** Assigning to a protected symbol is `Set::wrsym`
  ("Symbol Pi is Protected."); a rule for a protected tag is `Set::write`
  ("Tag Sin in Sin[x] is Protected."), as in the reference.
- **Rule order.** A rule goes ahead of the stored rules whose patterns are
  more general (see `DownValues`, Order). So `c[e_, {f___}, o:OptionsPattern[]]`
  is tried before `c[e_, f_, o:OptionsPattern[]]` whichever was defined first
  (v0.586; by pattern generality since v0.6022).
- **Orderless heads.** A DownValue whose first pattern slot names a head fires
  through any argument of an `Orderless` head: `Pair[n_ x_Momentum, y_] := ...`
  applies to `Pair[Momentum[p], lm Momentum[p]]` (v0.586).
- **Nested Part targets.** `x[[i]][[j]] = v` (and `+=`) is `x[[i, j]] = v` when
  the inner indices are integers (v0.586).
- **Redefinition.** A new rule replaces an existing one whose left-hand side
  (together with any `/;` condition) is the same up to renaming of pattern
  variables. The replacement covers both the pattern and the right-hand side:
  `f[x_] := x^2; f[y_] := y^3` leaves one rule, so `f[2]` is `8`. Rules with
  different conditions coexist (`f[x_] := a /; x > 0; f[y_] := b /; y < 0` keeps
  both). Each rule stores an alpha-invariant hash, so a definition costs one
  integer compare per stored rule unless the hashes agree (v0.475).
- A `PatternSequence`, `OrderlessPatternSequence`, `Longest`/`Shortest` of a
  sequence, or an `Alternatives` with a sequence branch in a definition's
  arguments makes its arity variable for DownValue dispatch:
  `q[PatternSequence[a_, b_], c_] := {a, b, c}` now fires for `q[1, 2, 3]`
  (v0.475).
- **Part assignment** (`x[[spec...]] = rhs`). A List `rhs` is distributed
  **level by level**: at a subscript that selects several positions (a list of
  positions, `All`, a `Span`), an `rhs` whose length equals the number of
  positions selected there is split across them, and each piece carries on to
  the deeper subscripts; any other `rhs` is stored whole at every position.
  A single-position subscript stores the whole `rhs` (v0.522):

  ```mathematica
  In[1]:= m = {{x, 3, 4}, {3, 4, 5}, {4, 5, 6}};
  In[2]:= m[[{1, 2}, {1, 2}]] = {{a, b}, {c, d}}; m
  Out[2]= {{a, b, 4}, {c, d, 5}, {4, 5, 6}}
  In[3]:= m[[All, 2]] = {p, q}; m
  Out[3]= {{a, {p, q}, 4}, {c, {p, q}, 5}, {4, {p, q}, 6}}
  ```

  A packed target is written on its buffer when the value fits (`m[[2]] =
  {0, 0, 0, 0}` keeps `m` packed), and otherwise through the List path,
  repacking when the result is again uniform. A visible `NDArray` coerces into
  its dtype; a value no buffer holds (a symbol, a ragged row) leaves it
  unchanged with `Set::ndarray`.

## AddTo (+=), SubtractFrom (-=), TimesBy (*=), DivideBy (/=)

Combine a value into an existing one and return the **new** value.

| operator | head | equivalent |
|---|---|---|
| `x += dx` | `AddTo[x, dx]` | `x = x + dx` |
| `x -= dx` | `SubtractFrom[x, dx]` | `x = x - dx` |
| `x *= dx` | `TimesBy[x, dx]` | `x = x dx` |
| `x /= dx` | `DivideBy[x, dx]` | `x = x/dx` |

All four have attribute `HoldFirst`. The target may be a symbol or a `Part`
expression referring to an existing value (`v[[2]] += 1`, `m[[i, j]] *= 2`). The
combination goes through `Plus`/`Times` and is then re-evaluated, so list
threading and symbolic simplification happen exactly as they would for the same
expression written out: `a = {1, 2}; a *= 3` gives `{3, 6}`, and `b = x; b *= 2`
gives `2 x`. If the target has no assigned value, an `X::rvalue` message is
emitted and the expression is left unevaluated.

```mathematica
In[1]:= q = 5; q *= 3; q
Out[1]= 15

In[2]:= w = {1., 2., 3.}; w[[3]] /= 4.; w
Out[2]= {1.0, 2.0, 0.75}
```

`Increment` (`++`), `Decrement` (`--`), `PreIncrement` and `PreDecrement` share
the same machinery; the first two return the **old** value, the rest the new one.

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
  (not its value) is operated on. `Protected` symbols are not affected.
- The *arguments* of `lhs` are evaluated exactly as `Set` evaluates them, so
  `i = 2; f[i] = 1; f[i] =.` removes the `f[2]` rule, and `sv[i][j] =.` finds a
  SubValue at every level. Patterns and Hold-attributed positions stay held.
- Always returns `Null`, whether or not a matching rule was found.

```mathematica
In[1]:= x = 5; x =.; x
Out[1]= x

In[2]:= f[x_] := x^2; f[x_] =.; f[3]
Out[2]= f[3]

In[3]:= fact[1] = 1; fact[n_] := n fact[n - 1]; fact[1] =.; DownValues[fact]
Out[3]= {HoldPattern[fact[n_]] :> n fact[n - 1]}
```

## UpSet (^=), UpSetDelayed (^:=)

> Gated behind the `UP_VALUES` build flag (on by default; `make UP_VALUES=0`
> compiles the whole subsystem out).

- `lhs ^= rhs` (`UpSet`) / `lhs ^:= rhs` (`UpSetDelayed`): assign `rhs` and
  associate the rule with **every distinct symbol at level one of `lhs`** — each
  argument that is a symbol, and the head of each argument. The rule becomes an
  **UpValue** of those symbols rather than a DownValue of the head of `lhs`.
- An UpValue of a symbol `s` fires when `s` appears at level one of an enclosing
  call — i.e. `s` is a direct argument, or the head of a direct argument — and is
  tried **before** that enclosing head's own DownValues or built-in.

**Features**:
- `UpSet` is `HoldFirst` (its `rhs` is evaluated once), `UpSetDelayed` is
  `HoldAll` (its `rhs` is held and re-evaluated on each use); both are `Protected`.
- The head of a pattern-constrained argument also counts: `a_mod + b_mod ^:= …`
  keys on `mod` (the head of the `Blank`), so the rule fires for `mod[…] + mod[…]`.
- An UpValue on a `Protected` symbol is refused with a message.

```mathematica
In[1]:= area[square] ^= s^2; UpValues[square]
Out[1]= {HoldPattern[area[square]] :> s^2}

In[2]:= prop[a, b[c]] ^= value; {UpValues[a], UpValues[b]}
Out[2]= {{HoldPattern[prop[a, b[c]]] :> value}, {HoldPattern[prop[a, b[c]]] :> value}}

In[3]:= f[g[x_]] ^:= h[x]; {f[g[2]], f[k[2]]}
Out[3]= {h[2], f[k[2]]}
```

## TagSet, TagSetDelayed, TagUnset (f /: …)

> Gated behind the `UP_VALUES` build flag.

- `f /: lhs = rhs` (`TagSet`) / `f /: lhs := rhs` (`TagSetDelayed`): assign `rhs`
  and associate the rule with the **single tag `f`**, which must occur in `lhs`.
  The rule is installed as an **own-, down- or up-value as appropriate**:
  - `f` *is* `lhs` → OwnValue (the tag is redundant: same as `f = rhs`).
  - `f` is the **head of `lhs`** → DownValue (same as `f[…] = rhs`).
  - `f` occurs as an **element of `lhs`, or the head of an element** (including
    the innermost head of a compound head, `f[x][y]`) → UpValue.
- `lhs` may be wrapped in `HoldPattern`, which is how a rule for an
  assignment head is usually written: `f /: HoldPattern[Set[f[a_], v_]] := ...`.
- `f /: lhs =.` (`TagUnset`): removes the matching own/down/up rule associated
  with the tag `f`. (Plain `lhs =.` does **not** reach an UpValue — the tag is
  required, matching the reference.)

**Features**:
- `TagSet`, `TagSetDelayed` and `TagUnset` are `HoldAll` and `Protected`; for the
  immediate `TagSet` the `rhs` is evaluated before the rule is stored.
- If the tag does not occur in `lhs`, `TagSet::tagnf` / `TagSetDelayed::tagnf`
  ("Tag k not found in lhs.") is issued and nothing is installed. If it occurs
  but too deep to key a rule -- `k /: h[sq[k[s_]]] := 1` -- the message is
  `::tagpos`: "Tag k in h[sq[k[s_]]] is too deep for an assigned rule to be
  found." On either refusal `TagSetDelayed` gives `$Failed` and `TagSet` its
  `rhs`. `SubValues` (`f[x][y]`, the tag as the head of the head) are not
  supported and are declined with a message.
- The rule remembers whether it was made by `=` or `:=`; `Definition` renders it
  back that way (see `Definition`).
- **Upvalues for the assignment heads.** An upvalue whose pattern is a `Set`,
  `SetDelayed` or `Unset` call intercepts that assignment. Nothing is stored;
  the rule's value is the result. Packages use this to keep bookkeeping when
  a value is assigned:
  ```mathematica
  sp /: Set[sp[a_, b_], z_] := (spv[a, b] = z);  sp[p, q] = m^2;  spv[p, q]  (* m^2 *)
  g /: Set[g, v_] := Print[v];  g = 3        (* prints 3; g gets no OwnValue *)
  ```
  A rule fires only when its whole pattern matches the assignment, so an
  upvalue with a different shape cannot intercept an unrelated `Set`. The
  definitional forms (`UpSet`, `TagSet`, ...) are not intercepted.
- A `MakeBoxes` or `Format` rule (`f /: MakeBoxes[f[x_], TraditionalForm] := ...`,
  `CLC /: MakeBoxes[CLC[x___][y___], fmt] := ...`) is stored in the tag's
  `FormatValues`, not its `UpValues`, as in the reference. See `FormatValues`.

```mathematica
In[1]:= mod /: mod[a_, p_] + mod[b_, p_] := mod[Mod[a + b, p], p];
        mod /: i_Integer mod[a_, p_] := mod[Mod[i a, p], p];
        mod[2, 5] + 3 mod[3, 5] - mod[1, 5]
Out[1]= mod[0, 5]

In[2]:= h /: f[h[x_]] := fh[x]; f[h[5]]
Out[2]= fh[5]

In[3]:= h /: f[h[x_]] =.; f[h[5]]
Out[3]= f[h[5]]
```

## Rules on Flat heads (Plus, Times, …)

A rewrite rule whose left-hand side has a **Flat** head — `Plus`, `Times`,
`Dot`, `**`, `And`, `Or`, or any symbol carrying the `Flat` attribute — matches a
**subset** of a longer operand sequence: the matched operands are rewritten and
the remaining ones are kept, re-wrapped in the head. This holds for UpValues,
DownValues, and `Replace`/`ReplaceAll` alike (matching the reference).

```mathematica
In[1]:= mod /: mod[a_, p_] + mod[b_, p_] := mod[Mod[a + b, p], p];
        mod[1, 5] + mod[2, 5] + mod[3, 5]
Out[1]= mod[1, 5]                 (* a pair is rewritten, then re-fires to a fixed point *)

In[2]:= mod[1, 5] + mod[2, 5] + x
Out[2]= mod[3, 5] + x             (* the non-matching term x is kept *)

In[3]:= 1 + x + y /. x + y -> z
Out[3]= 1 + z
```

- A match that **consumes every operand** wins over a subset ("leftover") match:
  with `pp /: Plus[pp, a_] := g[a]`, the input `pp + x + y` binds `a_` to the whole
  group `x + y` (→ `g[x + y]`), not `x` with `y` left over. The leftover match is
  taken only when no full-consumption match exists.
- `MatchQ` stays a **whole-expression** test and does *not* do leftover matching:
  `MatchQ[mod[1,5] + mod[2,5] + x, mod[a_,p_] + mod[b_,p_]]` is `False`, as in the reference.
- For a Flat head that is **not** `Orderless` (e.g. `Dot`, `**`, a user head with
  only `Flat`), the matched operands must be a **contiguous run in order**, but the
  run may sit anywhere — at the front, in the interior, or at the end. The
  **leftmost** run is taken, and a single application rewrites **one** run (as
  `ReplaceAll` rewrites each position once); `//.`, or a DownValue/UpValue
  re-firing to a fixed point, rewrites the rest. `Plus` and `Times` are
  `Orderless`, so there every subset of operands is reachable.
- **Cost.** Offsets and spans that cannot match are skipped before any work, using
  the pattern's first single-operand element (`op[j_]`, `Not[x_]`) as an anchor and
  rejecting a group that cannot equal an already-bound name (`ff[x_, x_]`) without
  building it. A failing `And[x_, Not[x_]]` over an n-term conjunction costs O(n)
  operand tests; `ff[x_, x_]` and `ff[x_, y_, x_]` cost about what the prefix-only
  match did before v0.462. The exception is a test on a group-absorbing first
  element, as in `ff[a_?Negative, b_]`: every (offset, span) group is a separate
  call to the test, so a failing rule of that shape makes about n²/2 test calls.
- `PatternSequence[...]` as an element of such a rule now also matches an interior
  run: `ff[x, a, b, y] /. ff[PatternSequence[a, b]] -> w` gives `ff[x, w, y]`.

```mathematica
In[1]:= x ** y ** z /. y ** z -> w
Out[1]= x ** w

In[2]:= SetAttributes[ff, Flat]; ff[x, y, z] /. ff[y, z] -> w
Out[2]= ff[x, w]

In[3]:= ff[a, b, a, b] /. ff[a, b] -> w
Out[3]= ff[w, a, b]               (* leftmost run, one rewrite *)

In[4]:= ff[a, b, a, b] //. ff[a, b] -> w
Out[4]= ff[w, w]

In[5]:= op /: op[i_] ** op[j_] := op[i + j];
        x ** op[1] ** op[2] ** y ** op[3] ** op[4]
Out[5]= x ** op[3] ** y ** op[7]
```

**v0.6032 (cost).** A failing rule on a wide Flat, Orderless sum no longer tries every subset of its operands. An Optional-coefficient element such as `c_. Cosh[a_]` takes only a single operand when its required part, a call on a fixed non-OneIdentity head, cannot equal a group. The leftover retry is skipped when the pattern already ends in a free `r___` that no condition reads. Matches are unchanged; 30-term sums are linear rather than exponential.

**v0.6121 (semantics and cost).** Under a Flat head, an element that is a call on another head takes exactly one operand, never a group of them, whether or not that head is OneIdentity and has defaulted arguments. This is the reference behaviour: `MatchQ[2 x y, y_^n_. x]` is `False` (it was `True`), and `ReplaceList[a + b + c, n_ o1_. + o2_ :> {n, o1, o2}]` binds `n` to `a`, `b` or `c` alone. The whole subject is still reached through the defaults: `ReplaceList[2 x, x_. y_^n_. :> {x, y, n}]` includes `{1, 2 x, 1}`. In an Orderless Flat pattern, a plain blank such as `o2_` or `o2_:0` is matched after the elements that take one operand, as `x__` already was. Together these make `n_?NumberQ o1_. + o2_:0` linear in the length of the sum; it had been about 3ⁿ (20 s at 14 terms). DownValue dispatch also keys every leading argument slot by its possible heads and call arity, and keys all slots of an Orderless head by assignment. A rule that cannot fit is skipped without running the matcher.

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
- Path components may be integers, the symbol `All` (selects every child at that level), or `Span` expressions such as `i ;; j` or `i ;; j ;; k` (including `UpTo[n]` bounds).
- Works on expressions with any head (not just `List`); after substitution the evaluator re-applies canonical ordering for `Orderless` heads such as `Plus` and `Times`.
- On an `Association` a position is a key, `Key[k]`, or a positional index over the entries, and the rules are tried against the *value*; `All`, `Span` and `0` work there too.
- The position list uses the same form as is returned by `Position`. `ReplaceAt[expr, rules, {}]` is an **empty list of positions** and replaces nothing; the position of the whole expression is the empty path, `{{}}`.
- A position that does not exist — an out-of-range index, an absent key, a malformed `Span`, or a path that runs into an atom — leaves `ReplaceAt` unevaluated, as `Part` does for `{a, b, c}[[99]]`.
- Position resolution is shared with [`MapAt`](functional-programming.md#mapat) via one walker (`expr_apply_at_path`, `src/part.c`), so the two agree on every position spec by construction.

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
- **A `Rational` or `Complex` number is an atom**: rules never reach its parts
  or its head (`{1/2} /. 2 -> 3` and `{1/2} /. Rational -> f` stay `{1/2}`, and
  `{I} /. Complex -> h` stays `{I}`). A rule for the number as a whole still
  matches: `{1/2} /. Rational[a_, b_] -> g[a, b]` is `{g[1, 2]}`. The same holds
  for `Replace` at any level and for `ReplaceRepeated`.
- **A guard on a delayed rule's right-hand side** is part of the match:
  `lhs :> body /; test` behaves as `lhs /; test :> body`, so a failing guard
  backtracks into other assignments of the pattern (which pair of a `Plus`)
  and then lets the next rule try; a scoped guard
  (`lhs :> Module[{t = ...}, v /; cond]`) is evaluated in its scope. The same
  holds for `Replace`, `ReplaceRepeated`, `ReplaceList` and `Cases`. A `Rule`
  (`->`) keeps a `Condition` value (v0.586).
- **Unevaluated replacements.** The outer `Unevaluated` of a replacement is
  removed as it is inserted: `x /. d -> Unevaluated[Sequence[]]` deletes `d`
  (also in `d^2`, giving `Power[2]` = `2`), and `Hold[d] /. d :> Unevaluated[1 + 1]`
  is `Hold[1 + 1]`. `ReplacePart` does the same (v0.586).
- Inside an association the rules reach the values only, never the keys:
  `<|a -> 1, b -> 2|> /. b -> a` is unchanged. An association is itself a rule
  set with literal keys: `{x, y} /. <|x -> 1|>` is `{1, y}` (see
  [Association atomicity](data-structures.md#association-atomicity-structural-functions)).

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

## ReplaceRepeated (//.)
`ReplaceRepeated[expr, rules]` (`expr //. rules`) applies `rules` repeatedly
with `ReplaceAll` semantics until the expression stops changing.
- A list of rule lists gives a list of results, one per rule list.
- **Option** `MaxIterations -> 65536`: the maximum number of scans, a
  non-negative integer or `Infinity`. It can be set with `SetOptions`.
- If the expression is still changing after `n` scans, the call emits
  `ReplaceRepeated::rrlim` ("Exiting after *expr* scanned *n* times", naming
  the original expression) and returns the expression reached. A fixed point
  reached on the `n`th scan still needs one more scan to confirm, so it
  messages too.
- An invalid `MaxIterations` value emits `ReplaceRepeated::iopnf`, and an
  unknown option emits `ReplaceRepeated::optx`. Either way the call stays
  unevaluated.
- `Protected`.

```
ReplaceRepeated[f[a], {f[a] -> f[b]}, MaxIterations -> 20]   (* f[b] *)
ReplaceRepeated[0, x_Integer :> x + 1, MaxIterations -> 5]   (* 5, with rrlim *)
```

## Dispatch
`Dispatch[{lhs1 -> rhs1, ...}]` is a rule set that `ReplaceAll` and
`ReplaceRepeated` index by left-hand side. Each subexpression then tries only
the rules that could match it, instead of every rule.
- `Replace`, `ReplaceList` and `ReplaceAt` accept it as their rules.
- A list of rule lists gives a list of results, as with plain rules.
- `Normal[Dispatch[rules]]` gives `rules`.
- `Dispatch[rule]` is `Dispatch[{rule}]`, and `Dispatch[Dispatch[r]]` is
  `Dispatch[r]`. Anything else stays unevaluated. `Protected`.

**How the index keys a rule:**
- A literal atom LHS (`a`, `1`, `"s"`) is keyed by that atom.
- A literal-headed LHS (`f[...]`) is keyed by its head.
- `HoldPattern`, `/;`, `?` and named-pattern wrappers are looked through, and
  `Verbatim[x]` is keyed by `x`.
- Every other rule is tried at every node: a pattern construct (`x_`, `_h`,
  `a | b`), a compound head, an association, or an LHS with an
  `Optional`/`x_.` argument (`x^n_.` matches a bare `x`).
- At each node the candidate rules are tried in their original order, so the
  rule that wins is always the one the plain list would pick.

**Performance:** for 5000 rules `vi -> i` over 20,000 symbols, the indexed
rule set takes 3.8 ms against 772 ms for the plain list.

**Divergence:** `Dispatch[...]` is a normal expression here, so `Length` and
`AtomQ` see its argument. The reference treats it as an atom, which prints as
`Dispatch[<n>]`.

## CompoundExpression (;)
- `expr1; expr2; ...`: Evaluates a sequence of expressions, returning the last one.

## DownValues, OwnValues

- `DownValues[f]`: the list of transformation rules associated with the symbol
  `f` as a function head — the rules set by `f[patt] = rhs` / `f[patt] := rhs`.
- `OwnValues[f]`: the list of rules for `f` used as a bare symbol — the value set
  by `f = rhs` / `f := rhs`.

**Features**:
- Each rule is returned as `HoldPattern[lhs] :> rhs` (a `RuleDelayed` whose
  left-hand side is wrapped in `HoldPattern`), matching the reference. The
  `HoldPattern` keeps the stored pattern from re-matching its own definition, and
  the `RuleDelayed` keeps the right-hand side unevaluated, so the returned list is
  **inert**: querying `DownValues[f]` for a recursive `f` does not re-trigger the
  recursion, and a right-hand side such as `RandomInteger[...]` is shown held
  rather than evaluated afresh on each query.
- Returns `{}` for a symbol that carries no rules of the requested kind.
- **Order** (v0.6022, as in the reference). The list is ordered by pattern
  generality, not by a score.
  - **Pattern-free definitions** (`f[1]`, `f[g[a]]`) come first of all, in
    definition order. They are also kept in a per-symbol hash index (v0.6023),
    so looking one up or redefining it does not scan the other stored values.
    A memoised `f[n_] := f[n] = ...` over 20000 values takes ~0.03 s; it took
    ~38 s, because every call matched against every stored value. A `Flat` or
    `Orderless` head still scans, since a literal `f[a, b]` matches more than
    its equal there.
  - **Strictly more general.** A new pattern rule goes ahead of the first
    stored rule that matches everything it matches, and more. So `f[x_]`
    follows `f[x_Integer]`, `f[x_?EvenQ]` and `f[x_ /; x > 0]`, and
    `f[x__]` follows `f[x_, y_]`.
  - **Incomparable.** Two patterns, neither of which covers the other, keep
    definition order however specific either looks:
    `f[x_Integer]; f[x_?EvenQ]` and `f[h[x__, y_]]; f[h[x_, y__]]`. So does
    FeynCalc's `h[b___G, di_G, c__G] /; ..` before a later
    `h[b___, G[5], c:G[_[__]].., d___]`, which a score had once moved ahead.
  - **Mutual cover.** Of two patterns that cover each other the newer goes
    first: `x_` and `x_^n_.`, or `h[x___, y___]` and `h[z___]`. Two that are
    the same up to renaming are the exception: conditional definitions
    `f[x_] := a /; c` kept side by side stay in definition order.
  - **Guards** narrow a pattern but do not order it otherwise:
    `f[p[a_, a_]] := 1; f[p[a1_, b1_]] := 2 /; True` tries the first rule
    first, and `f[x_] := 1; f[x_ /; x > 0] := 2` gives `f[3]` = `2`.
  - The comparison is structural and conservative. It understands typed
    blanks, sequences, `Repeated` (`x:p..` is narrower than `x__`), defaults
    (`x_^n_.` also matches what `x_` matches), alternatives, tests and
    nonlinear patterns (`f[x_, x_]` is never taken as the more general).
    When it cannot tell, it keeps definition order.

```mathematica
In[1]:= f[x_] := x^2; DownValues[f]
Out[1]= {HoldPattern[f[x_]] :> x^2}

In[2]:= a = 5; OwnValues[a]
Out[2]= {HoldPattern[a] :> 5}
```

## Format, FormatValues, MakeBoxes

- `Format[expr] := rhs` / `Format[expr, form] := rhs`: display `expr` as `rhs`.
- `FormatValues[f]`: the display rules of `f`, as `HoldPattern[lhs] :> rhs`.
- `FormatValues[f] = {rules}` replaces them; `FormatValues[f] =.`, `Clear[f]`
  and `ClearAll[f]` remove them.
- `f /: MakeBoxes[lhs, form] := boxes` and `f /: Format[lhs] := rhs` store
  their rule in `FormatValues[f]`, never in `UpValues[f]`.

**Features**:
- An untagged `Format` is keyed on the symbol of `expr`, or on its innermost
  head. The stored rules are the reference's:
  `Format[k[x_]] := x + 1` gives `{HoldPattern[MakeBoxes[k[x_], FormatType_]] :>
  Format[x + 1, FormatType], HoldPattern[k[x_]] :> x + 1}`.
  `Format[k[x_], OutputForm] := ...` stores only the plain rule, and other
  forms store only the `MakeBoxes` rule.
- `ToString` (OutputForm) and the REPL's `Out[]` apply the plain rules, inside
  `Hold`, `HoldForm` and `HoldPattern` too. A result is formatted again
  (`Format[a[x_]] := b[x]; Format[b[x_]] := 9` shows `a[1]` as `9`). A rule
  that keeps producing a match stops with `Format::forml`.
- `InputForm`, `FullForm` and evaluation never see the rules. A script's
  `Print` writes the raw expression, as `wolframscript` does.
- Mathilda has no box typesetting, so `MakeBoxes` rules are stored and listed
  but not rendered.
- `Format` is `Protected`, `FormatValues` is `HoldAll`/`Protected`, and
  `MakeBoxes` is `HoldAllComplete` (a `System` symbol). An untagged
  `MakeBoxes[...] := ...` is a DownValue of `MakeBoxes`. Assigning a format to
  a `Protected` symbol gives `SetDelayed::write`.

```mathematica
In[1]:= Format[k[x_]] := x + 1; {ToString[k[3]], ToString[k[3], InputForm], k[3] === k[3]}
Out[1]= {4, k[3], True}

In[2]:= f /: MakeBoxes[f, TraditionalForm] := "F"; {UpValues[f], FormatValues[f]}
Out[2]= {{}, {HoldPattern[MakeBoxes[f, TraditionalForm]] :> F}}
```

## ClearAll, Remove

- `ClearAll[s1, s2, ...]`: Clears all values (OwnValues and DownValues),
  attributes and usage messages for the named symbols, leaving the symbols
  themselves in the symbol table. `ClearAll[{s1, s2, ...}]` accepts a list of
  specs.
- `Remove[s1, s2, ...]`: Removes the named symbols completely, deleting their
  definitions from the symbol table. A later reference recreates a fresh, empty
  symbol of the same name. `Remove[{s1, s2, ...}]` accepts a list of specs.

**Features**:
- `ClearAll` and `Remove` both have attributes `{HoldAll, Protected}`. Both hold
  their arguments, so they operate on the symbol, not its current value.
- Neither affects symbols with the attribute `Protected`. This is what prevents
  `Remove`/`ClearAll` from ever deleting or wiping a built-in.
- `ClearAll`, unlike `Clear`, also removes attributes and the usage message.
- Both return `Null`.
- **String specs are name patterns.** A string argument is resolved the way
  `Names` resolves it: `*` matches any characters, `@` matches characters that
  are not uppercase letters, and a context prefix restricts the match. So
  `Remove["Global`x*"]`, `Remove["ctx`*"]`, `ClearAll["b*"]` and
  `Clear["c*"]` act on every symbol `Names[spec]` lists. A plain name, with or
  without its context, names that one symbol.
- `Remove` emits `Remove::rmnsm` when a string spec matches nothing, and
  `Remove::rmptc` for each `Protected` symbol it skips.

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
- The typical sequence for adding rules to an existing symbol is
  `Unprotect[f]; definition; Protect[f]`.
- This applies to builtins too: `Unprotect[Sin]` really removes `Protected`, so
  `Attributes[Sin]` is `{Listable, NumericFunction}` and `Sin[x] = ...` is
  accepted until `Protect[Sin]`. `SetAttributes`/`ClearAttributes` change a
  builtin's attributes in the same way (issue #103).

```mathematica
In[1]:= f[x_] := x^2; Protect[f]
Out[1]= {"f"}

In[2]:= f[x_, y_] := x + y
(* SetDelayed::wrsym: Symbol f is Protected. *)

In[3]:= Unprotect[f]
Out[3]= {"f"}

In[4]:= Protect[f]; Protect[f]
Out[4]= {}

In[5]:= Check[f[x_] := 0, "refused"]
Out[5]= "refused"
```

## Options, SetOptions, OptionValue

Default option settings are stored per symbol (the `DefaultValues`-equivalent),
separate from its OwnValues/DownValues. Built-in option-accepting functions ship
with their defaults pre-registered, so `Options[f]` is queryable for both kernel
and user symbols.

- `Options[f]` (symbol): the list of default option rules `{name -> value, ...}`,
  or `{}` if none.
- `Options[expr]` (compound expression): the option rules explicitly present in
  the expression's arguments — e.g. the options of a `Graphics` object.
- `Options[obj, name]` / `Options[obj, {n1, n2, ...}]`: the selected setting(s),
  as a list of rules.
- `Options[f] = {name -> value, ...}`: redefines **all** of `f`'s default options
  at once (works through `Set`, bypassing the `Protected` attribute on `Options`).
- `SetOptions[s, name -> value, ...]`: changes individual existing defaults and
  returns the new `Options[s]`. It works on `Protected` symbols. It cannot
  **add** an option — an unknown name raises
  `SetOptions::optnf` and leaves the call unevaluated. Use
  `AppendTo[Options[s], name -> value]` (or `PrependTo`) to add one.
- `SetOptions[s, {name -> value, ...}]` takes the rules as a list, and
  `SetOptions[{s1, s2, ...}, rules]` sets them on every symbol in the list and
  returns the list of their new options:
  `SetOptions[{g, h}, opt -> False]` is `{{opt -> False, o2 -> 1}, {opt -> False}}`.
- `OptionValue[name]`, `OptionValue[f, name]`, `OptionValue[f, opts, name]`,
  `OptionValue[f, opts, name, Hold]`: resolve a single option value from the
  explicit rules `opts` first, then the defaults derived from `f` (a symbol, a
  rule, or a list of either). The bare and two-argument forms resolve inside a
  rule whose left-hand side used `OptionsPattern[]` (see
  [Pattern Matching](pattern-matching.md)). Option names may be symbols or
  strings interchangeably, ignoring contexts; a trailing `Hold` wraps the
  resolved value in `Hold`. Lists of rules may nest (`OptionValue[f, {{a -> 1}},
  a]` is `1`), searched in order (v0.592).

**Features**:
- `Options`, `SetOptions`, and `OptionValue` all have attribute `{Protected}`.
- Default options survive `Clear[f]` (only rules are cleared) and are removed
  with the symbol by `Remove[f]`.
- The registered defaults mirror the options each builtin's evaluator actually
  honors, with the value it falls back to when the option is absent. The sweep
  covers, among others: `Integrate` (`Method`), `Limit` (`Direction`,
  `Assumptions`), `Series`/`PowerExpand` (`Assumptions`), `D`/`Dt`
  (`NonConstants`), `Sum`/`Product` (`Method`), `Simplify`/`FullSimplify`
  (`Assumptions`, `ComplexityFunction`, `TransformationFunctions`),
  `GroebnerBasis` (`MonomialOrder -> Lexicographic`, `CoefficientDomain ->
  Rationals`, `Method`, `Sort`, `Modulus`), `Factor`/`Together`/`Cancel`/`Apart`
  and the `Polynomial*` family (`Extension`, `Modulus`), `IrreduciblePolynomialQ`
  /`SquareFreeQ` (`GaussianIntegers`, ...), the eigen/linear-algebra heads
  (`Eigenvalues`/`Eigenvectors` → `Cubics`, `Quartics`; `LinearSolve`,
  `LeastSquares`, `NullSpace`, `MatrixRank`, `PseudoInverse`,
  `SingularValueDecomposition`), `PrimeQ`/`CoprimeQ`/`FactorInteger`
  (`GaussianIntegers`), the `Random*` heads (`WorkingPrecision`), the numerical
  calculus/root-finding heads (`NIntegrate`, `NSum`, `NProduct`, `NLimit`, `ND`,
  `NSeries`, `NResidue`, `FindRoot`, `FindMinimum`, `FindMaximum`, `NRoots`,
  `NSolve`, `Solve`, `Fit`), `Plot`, and the structural family that reads
  `Heads` (`Cases`/`Count`/`DeleteCases`/`MemberQ`/`Map`/`Apply`/`MapAll`/`Level`
  /`Depth`/`LeafCount` default `Heads -> False`; `Position`/`FreeQ` default
  `Heads -> True`).

```mathematica
In[1]:= Options[LinearSolve]
Out[1]= {Method -> Automatic, Modulus -> 0, ZeroTest -> Automatic}

In[2]:= Options[f] = {a -> 1, b -> 2}; f[OptionsPattern[]] := {OptionValue[a], OptionValue[b]}

In[3]:= {f[], f[a -> 17], f[b -> 18], f[a -> 17, b -> 18]}
Out[3]= {{1, 2}, {17, 2}, {1, 18}, {17, 18}}

In[4]:= SetOptions[f, b -> 20]; f[]
Out[4]= {1, 20}

In[5]:= SetOptions[f, c -> 3]
(* SetOptions::optnf: c is not a known option for f. *)
Out[5]= SetOptions[f, c -> 3]

In[6]:= AppendTo[Options[f], c -> 3]
Out[6]= {a -> 1, b -> 20, c -> 3}
```

**`OptionValue` in a condition.** A rule whose LHS carries `OptionsPattern[]`
resolves `OptionValue[name]` in its `/;` guard exactly as it does in its RHS,
against the supplied options and then the head's defaults:
```
Options[g] = {o -> 1};
g[x_, OptionsPattern[]] := x /; OptionValue[o] == 1
{g[3], g[3, o -> 2]}                       (* {3, g[3, o -> 2]} *)
```
A failing guard falls through to the next rule. This also applies to a guard
written on the LHS (`g[x_, OptionsPattern[]] /; OptionValue[o] == 1 := ...`).
That case is deliberately more permissive than the reference, which leaves
such a call unevaluated.

- **An aliased option name** resolves at lookup: after `Options[f] = {al -> 1}`
  and a later `al = real`, `OptionValue[al]` (which arrives as `real`) still
  finds the option (v0.586).

## OptionQ
`OptionQ[e]` gives `True` if `e` is an option rule, or a possibly nested list of
them, and `False` otherwise.
- An option rule is a `Rule` or `RuleDelayed` whose name is a symbol or a
  string: `a -> 1`, `"a" :> 1`. A rule with any other name (`1 -> 2`,
  `f[x] -> 2`) is not one.
- Lists may nest and may be empty: `{}`, `{{a -> 1}, {b -> 2}}`.
- An `Association` is `False`.
- `Protected`. A wrong argument count emits `OptionQ::argx`.

The usual use is as a pattern test that collects trailing options:
```
f[x_, opts___?OptionQ] := {x, opts}
f[1, a -> 2]                      (* {1, a -> 2} *)
```

## FilterRules
`FilterRules[rules, spec]` gives the rules whose left-hand sides match `spec`,
in their original order. Duplicates are kept.

| `spec` | Selects |
|---|---|
| a symbol or string | rules with that name, ignoring context; `"b"` and `b` select the same rule |
| `{s1, s2, ...}` | rules any `si` selects |
| a rule, or a list of rules such as `Options[f]` | rules with those rules' names |
| `Except[s]`, `Except[c, p]` | the complement |
| any other pattern | rules whose name matches it, e.g. `_`, `a \| b`, `x_ /; cond` |

- `rules` may be a single rule or nested lists of rules.
- A non-rule element emits `FilterRules::rep` and leaves the call unevaluated.
- `Protected`.

```
FilterRules[{a -> 1, b -> 2, c :> 3}, {a, c}]     (* {a -> 1, c :> 3} *)
FilterRules[{opts}, Options[Plot]]                 (* the options Plot accepts *)
FilterRules[{opts}, Except[Options[Plot]]]         (* the rest *)
```
