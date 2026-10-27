# Pattern Matching

- `_`: `Blank[]`.
- `__`: `BlankSequence[]`.
- `___`: `BlankNullSequence[]`.
- `x_`: `Pattern[x, _]`.
- `_h`: Matches expression with head `h`. Head-typed blanks cover the atomic
  numeric heads too: `_Rational` matches `1/2` (`Rational[1, 2]`), `_Complex`
  matches `3 + 4 I` (`Complex[3, 4]`), and `_Integer` matches an arbitrary-size
  bigint.
- `p ? test`: `PatternTest[p, test]`. On a sequence blank (`__?test`,
  `___?test`) the predicate is applied to every matched element.
- `lhs /; condition`: Conditional matching. Chained `/;` nests right, i.e.
  `p /; c1 /; c2` is `Condition[p, Condition[c1, c2]]`, and every leaf condition
  must hold (`Condition` is `HoldAll`, so the leaves are evaluated individually).
  A failing condition (or a failing later argument) **backtracks into nested
  patterns** — every assignment of a nested `Orderless` / `Flat` / sequence match is
  tried, as in the reference (v0.451): `MatchQ[g[z x^k], g[a_ f_] /; FreeQ[a, x]]` is
  `True` (`a -> z`), where the matcher previously committed to the nested match's
  first split (`a -> x^k`) and failed.
- `p1 | p2 | ...`: `Alternatives`, matches any one of the `pi`. A name that
  appears only in a branch that was not taken stays unbound during the match,
  and is bound to `Sequence[]` once the whole match succeeds, so it stands for
  nothing on the right-hand side: `1 /. (x_ | y_) :> {x, y}` is `{1}`. A later
  occurrence can still bind the name: `f[1, 2] /. f[x_ | y_, y_] :> {x, y}` is
  `{1, 2}` (v0.475).
- `sym:p`: `Pattern[sym, p]`. The `:` binds only the token just before it (see
  [Operators](../operators.md)), so `n_.*a:_g` is `n_. * (a:_g)`. The name is
  bound first and `p` is then matched with the rest of the match as its
  continuation, so a later failure backtracks into `p`'s other assignments:
  `ReplaceList[2 x, p:x_.*y_ :> {x, y}]` lists all three matches (v0.475).
- `Except[c]`, `Except[c, p]`: matches anything not matching `c` (and, in the
  two-argument form, also matching `p`).
- `HoldPattern[p]`: matches as `p` but keeps `p` unevaluated when the rule is
  formed (so e.g. `x_ + x_` is preserved as a sum pattern rather than
  collapsing to `2 x_`).
- `Verbatim[expr]`: matches `expr` taken literally — the pattern constructs
  inside `expr` are not interpreted, so `Verbatim[x_]` matches only the literal
  `x_` (`Pattern[x, Blank[]]`).
- `PatternSequence[p1, ..., pm]`: matches a run of arguments each matching
  `p1..pm`; `PatternSequence[]` matches zero arguments, and `x:PatternSequence[
  ...]` binds `x` to the matched `Sequence`.
- `OrderlessPatternSequence[p1, ..., pm]`: matches a **contiguous run** of
  arguments at its own position, as any sequence element does, with the `pi`
  matching the run's elements in any order (even under a non-`Orderless` head).
  `f[2, 1, b]` matches `f[OrderlessPatternSequence[1, 2], ___]` but `f[1, b, 2]`
  does not. A `pi` may be a sequence pattern (`f[OrderlessPatternSequence[x__,
  c]]` matches `f[a, b, c]` with `x -> Sequence[a, b]`). Literal `pi` are placed
  first, on the first equal element. `x:OrderlessPatternSequence[...]` binds `x`
  to the run in its original order. Before v0.475 the `pi` could take any
  argument anywhere, which matched too much (`f[1, 2, 3]` matched
  `f[3, OrderlessPatternSequence[1, 2]]`).
- `Longest[p]` / `Shortest[p]`: prefer the longest / shortest match for a
  sequence pattern `p`; transparent at a non-sequence position.
- A sequence pattern used at the **top level** matches the single subject as a
  length-1 sequence: `MatchQ[a, a..]`, `MatchQ[a, (_)..]` and `MatchQ[a, a...]`
  are all `True` (as are `__`/`___` there). A `List` is one expression, so it
  is not decomposed against a bare top-level sequence pattern —
  `MatchQ[{a, b}, (a | b)..]` is `False` (use `{(a | b)..}`). Nested repeats
  such as `Repeated[RepeatedNull[a]]` (`(a...)..`) are matched by routing each
  repetition through this top-level rule.

The matcher's recursion depth is bounded by the same user-visible
`$RecursionLimit` the evaluator uses (default `1024`): a subject or pattern
nested deeper than that yields a graceful non-match rather than overflowing the
C stack. Raise `$RecursionLimit` to match deeper structures.

## Cases
Gives a list of elements that match a pattern.
- `Cases[{e1, e2, ...}, pattern]`
- `Cases[{e1, ...}, pattern -> rhs]`
- `Cases[expr, pattern, levelspec]`
- `Cases[expr, pattern -> rhs, levelspec]`
- `Cases[expr, pattern, levelspec, n]`
- `Cases[pattern]`

**Features**:
- Uses standard level specifications. A positive bound counts levels from the
  top and a negative bound `-k` counts by `Depth`. The two combine
  independently, so a mixed spec such as `{2, -1}` (level 2 down to the leaves)
  works in `Cases`, `Count`, `Position`, `DeleteCases`, `MemberQ` and `FreeQ`.
  Before v0.475 it matched nothing.
- Defaults to level `{1}`.
- A plain `_h` pattern is tested by comparing heads, skipping the matcher, on
  every compound node of a level traversal (v0.6052).
- Option `Heads -> True` includes heads of expressions.
- Operates in depth-first post-order traversal.
- `Cases[expr, p :> v /; test]` takes the guard as part of the match, as
  `ReplaceAll` does: `Cases[{1, 2, 3}, x_ :> (x^2 /; x > 1)]` is `{4, 9}` (v0.586).

## DeleteCases
Removes the elements of an expression that match a pattern.
- `DeleteCases[expr, pattern]`
- `DeleteCases[expr, pattern, levelspec]`
- `DeleteCases[expr, pattern, levelspec, n]`
- `DeleteCases[pattern]`

**Features**:
- Uses standard level specifications, defaulting to level `{1}`.
- Option `Heads -> True` tests heads as well; deleting a head is equivalent to applying `FlattenAt` at that point, splicing the remaining arguments into the parent.
- Traverses `expr` in depth-first post-order (leaves before roots) so that the `n` budget is spent on deeper matches before shallower ones.
- The match test is applied to the original subexpression (not the version with children already deleted), matching the reference semantics.

## Position
Gives a list of the positions at which objects matching a pattern appear.
- `Position[expr, pattern]`
- `Position[expr, pattern, levelspec]`
- `Position[expr, pattern, levelspec, n]`
- `Position[pattern]`

**Features**:
- Defaults to levels `{0, Infinity}` with `Heads -> True`.
- Yields lists of indices in lexicographic order.

## FirstPosition
Gives the position of the first element matching a pattern, in depth-first order.
- `FirstPosition[expr, pattern]` — the position, or `Missing["NotFound"]` if none matches.
- `FirstPosition[expr, pattern, default]` — `default` if no match is found.
- `FirstPosition[expr, pattern, default, levelspec]`

**Features**:
- Attributes `{HoldRest, Protected}`; `default` is evaluated only when it is returned.
- Same defaults as `Position`: levels `{0, Infinity}` with `Heads -> True`; a position of `{}` is the whole of `expr`.
- Over an association the position is a key, e.g. `{Key[k], ...}`.
- Delegates to `Position` (with the first-match cap), so traversal, level specs, the `Heads` option, and association handling match `Position` exactly.

```mathematica
In[1]:= FirstPosition[{a, b, c, b}, b]
Out[1]= {2}

In[2]:= FirstPosition[{{1, 2}, {3, 4}}, 4]
Out[2]= {2, 2}

In[3]:= FirstPosition[{1, 2, 3}, 5, None]
Out[3]= None
```

## Count
Gives the number of elements or subexpressions that match a pattern.
- `Count[list, pattern]`
- `Count[expr, pattern, levelspec]`
- `Count[pattern]`

**Features**:
- Uses standard level specifications, defaulting to level `{1}`.
- Option `Heads -> True` evaluates heads of expressions.
- Effectively returns the length of the list generated by `Cases` for the same arguments.
- `p:def`: `Optional[p, def]`, matches `p`, but if omitted defaults to `def`.
- `x_.`: `Optional[Pattern[x, _]]`, matches `x_`, but if omitted defaults to globally specified `Default` values.
- `p..`: `Repeated[p]`, matches a sequence of one or more expressions, each matching `p`.
- `p...`: `RepeatedNull[p]`, matches a sequence of zero or more expressions, each matching `p`.
- `Repeated[p, max]`, `Repeated[p, {min, max}]`, `Repeated[p, {n}]`: match sequences of length bounded by `max` or `{min, max}` (`max` may be `Infinity`). Equivalent properties exist for `RepeatedNull`. The bounds count *repetitions*, so they are argument counts only when `p` matches one argument at a time. A malformed spec (negative or non-integer count, `{min, max}` with `min > max`, `Repeated[p, 0]`) issues `Repeated::range` or `Repeated::order` and matches nothing. A bound too large for any argument list (`10^30`) acts as unbounded. `Repeated` and `RepeatedNull` are `Protected`.
- **Tokenization.** `..` and `...` are single tokens that take precedence over a decimal point: `0..` is `Repeated[0]`, `0...` is `RepeatedNull[0]`, and `1.5..` is `Repeated[1.5]`. After a headed or sequence blank the dots are the repeat operator (`_Integer..`, `x__..`). After a bare single blank `_.` is formed first, as in the reference, so `_...` is `(_.)..` and `_..` is a syntax error; write `(_)..` (that is also how it prints).
- **Unnamed sub-patterns** inside a repeat must bind the same value in every repetition: `{1, 1}` matches `{(x_ /; x > 0)..}`, `{1, 2}` does not. An unconstrained `x:p..` followed by more patterns takes the *shortest* run (`{1, 1, 2} /. {x:1.., y___} :> {{x}, {y}}` is `{{1}, {1, 2}}`); wrap it in `Longest` for the longest run.
- **Backtracking into units.** When the rest of the match fails, or
  `ReplaceList` asks for more matches, the unit search tries the next way to
  split the run: `ReplaceList[f[a, a], f[x__..] :> {x}]` is `{{a}, {a, a}}`
  (v0.475).
- **Multi-argument units.** The repeated pattern may span several arguments: `{a, b, a, b}` matches `{PatternSequence[a, b]..}`, `{c, a, b}` matches `{(PatternSequence[a, b] | c)..}`, and `{1, 2, 1, 2} /. {x__..} :> {x}` is `{1, 2}`. A unit of fixed width is checked in one linear pass. A variable-width unit is searched over its feasible widths without recursion, so a 40 000-argument subject stays linear even when the match fails.
- **Sequence alternatives.** An `Alternatives` whose branches span other than one argument matches as a sequence, named or not: `{a, b, d}` matches `{PatternSequence[a, b] | c, d}`, and `{a, b, d} /. {x:(PatternSequence[a, b] | c), y___} :> {{x}, {y}}` is `{{a, b}, {d}}`.

## Default Values
- `Default[f] = val`: Sets the default value for arguments of `f`.
- Built-in default values include `0` for `Plus` and `1` for `Times` and `Power`.
- Used automatically by the pattern matcher when `_.` is present in the pattern.

## Attributes in matching
- **Orderless.** Pattern arguments may take the subject's arguments in any
  order. A literal argument (no pattern construct) is matched to the first
  equal argument only, because every equal choice leaves the same bindings:
  `ReplaceList[o[a, a], o[a, x_] :> x]` is `{a}`.
- **OneIdentity.** `f[p1, ..., pn]` matches a lone `x` as though it were
  `f[x]` only through **defaults**: some `pi` must be Optional, so that the
  others match `x` while it takes its default (`x` matches `Times[n_., x_]`
  with `n -> 1`). Without an Optional, `f[x_]` matches neither `x` nor
  `f[a, b]`, and `f[]` under Flat+OneIdentity does not match `f[x_]` (before
  v0.475 all three matched). Under Flat, a subject that already has head `f`
  is covered by the ordinary path, so its default match is listed once.
  A guarded form keeps this: `MatchQ[h[g], a_. h[b_] /; FreeQ[{a, b}, h]]` is
  `True` with `a -> 1` (before v0.586 the guard made a head mismatch final).
- **A call with a pattern head matches one operand of a Flat call**
  (v0.6020), as in the reference. `g_[y__]` and `_[x_]` never take a group
  `h[...]` of several operands: `ReplaceList[ff[a[1], b[2], c[3]],
  ff[g_[y__], c[3]] :> g]` is `{}` for a Flat, Orderless `ff`, and
  `ff[a[1], b[2], c[3]] /. ff[g_[y__], c[3]] :> {g, {y}}` is
  `ff[b[2], {a, {1}}]`. Before, `g` could bind `ff` over `a[1], b[2]`. In a
  product, the UpValue `q /: q[z_, x_] f_[y__] := ...` used to try
  `f -> Times` over all the other factors.
- **Flat without OneIdentity** differs from the reference in `ReplaceList` enumeration.
  The reference binds a one-operand slice both as `h[b]` and as `b` and repeats matches;
  Mathilda binds `b` once. These cases are recorded, report-only, in
  `tests/match_exhaustive_observe.m`.

## OptionsPattern
- `OptionsPattern[]` matches the trailing option arguments of a call — each a
  `name -> value` (or `name :> value`) rule, or a `List` of such rules, which is
  flattened. As the least-specific, variable-arity slot it greedily consumes the
  whole trailing option run; any remaining non-option argument makes the rule
  fail to match.
- Inside the rule's right-hand side, `OptionValue[name]` (and
  `OptionValue[f, name]`) read the matched options, falling back to the defaults
  registered with the head via `Options` (see
  [Assignment and Rules](assignment-and-rules.md)).
- `OptionsPattern[{name -> default, ...}]` (or a single rule) carries its own
  defaults, which `OptionValue` uses for an option that was not passed:
  `q[x_, OptionsPattern[{o -> 1}]] := OptionValue[o]` gives `q[1] == 1` and
  `q[1, o -> 2] == 2` (v0.475).
- Option lists may nest to any depth: `gg[x_, OptionsPattern[]]` matches
  `gg[1, {{a -> 5}}]` (v0.592).
- **Any position** (v0.6011). `OptionsPattern[]` need not be the last element.
  Before later elements it takes a leading run of option arguments, longest
  first and down to none, and the rest go to the later elements:
  `MatchQ[t[x], t[OptionsPattern[], a_]]` is `True`, and
  `t[a -> 1, x] /. t[o:OptionsPattern[], b_] :> {{o}, b}` is `{{a -> 1}, x}`.
  On an `Orderless` head the slot sorts first: `s[a__, OptionsPattern[]]`
  becomes `s[OptionsPattern[], a__]` and still matches `s[x]`.
- **Named and guarded.** `opts:OptionsPattern[]` binds `opts` to the option
  arguments as written: a list of rules stays a list (`f[x_, o:OptionsPattern[]]
  := {o}; f[1, {a -> 1}, b -> 2]` is `{{a -> 1}, b -> 2}`). A guard on it,
  `(opts:OptionsPattern[]) /; test`, is checked once the options are bound
  (`Isolate[ex_, (opts:OptionsPattern[]) /; opts =!= {}] := ...`); it used to
  never match (v0.591).

```mathematica
In[1]:= Options[g] = {a -> 1, b -> 2};
        g[OptionsPattern[]] := {OptionValue[a], OptionValue[b]}
        {g[], g[a -> 17], g[{a -> 3, b -> 4}]}
Out[1]= {{1, 2}, {17, 2}, {3, 4}}
```

