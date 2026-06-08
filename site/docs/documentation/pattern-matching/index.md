# Pattern Matching

15 built-in function(s) in this category.

- [`Blank`](Blank.md) — _ or Blank[] represents any single expression.  _(Stable)_
- [`BlankNullSequence`](BlankNullSequence.md) — ___ or BlankNullSequence[] represents a sequence of zero or more expressions.  _(Stable)_
- [`BlankSequence`](BlankSequence.md) — __ or BlankSequence[] represents a sequence of one or more expressions.  _(Experimental)_
- [`Cases`](Cases.md) — Cases[{e1, e2, ...}, pattern] gives a list of the ei that match the pattern.  _(Stable)_
- [`Count`](Count.md) — Count[list, pattern] gives the number of elements in list that match pattern.  _(Stable)_
- [`Default`](Default.md) — Default[f]  _(Experimental)_
- [`DeleteCases`](DeleteCases.md) — DeleteCases[expr, pattern] removes all elements of expr that match pattern.  _(Stable)_
- [`HoldPattern`](HoldPattern.md) — HoldPattern[expr]  _(Stable)_
- [`Longest`](Longest.md) — Longest[p] is a pattern object that matches the longest sequence consistent with the pattern p.  _(Stable)_
- [`MatchQ`](MatchQ.md) — MatchQ[expr, form]  _(Stable)_
- [`Optional`](Optional.md) — patt:def or Optional[patt, def]  _(Stable)_
- [`Position`](Position.md) — Position[expr, pattern] gives a list of the positions at which objects matching pattern appear in expr.  _(Stable)_
- [`Repeated`](Repeated.md) — p.. or Repeated[p] is a pattern object that represents a sequence of one or more expressions, each matching p.  _(Stable)_
- [`RepeatedNull`](RepeatedNull.md) — p... or RepeatedNull[p] is a pattern object that represents a sequence of zero or more expressions, each matching p.  _(Stable)_
- [`Shortest`](Shortest.md) — Shortest[p] is a pattern object that matches the shortest sequence consistent with the pattern p.  _(Stable)_
