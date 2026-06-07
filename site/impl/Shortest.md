---
source: src/match.c
---
`Shortest[p]` is a pattern-matching modifier resolved inside the argument matcher in `src/match.c`, with no builtin. When peeling a leading pattern, the matcher unwraps wrapper heads (`Pattern`, `Optional`, `Longest`, `Shortest`) and sets an `is_shortest` flag (a following `Longest` overrides it). The flag biases the backtracking sequence matcher to prefer the *fewest* arguments consistent with the rest of the pattern (e.g. for `Optional`+`Shortest` it tries the absent/default binding first), inverting the default greedy ("longest") preference. It does not change which matches are *possible*, only the order they are tried.
