---
source: src/eval.c
---
**Algorithm.** `Nothing` is not a builtin function but a special symbol handled
by the evaluator: it is the identity element for list construction and vanishes
from any `List` it appears in. `strip_nothing` (`src/eval.c`), called during
`evaluate_step`, scans a `List`'s arguments and removes every element that is
either the symbol `Nothing` or any `Nothing[...]` form (an `EXPR_FUNCTION` whose
head is `Nothing`). It makes one counting pass to find how many survive, and only
if some element is a `Nothing` does it rebuild the argument array — freeing the
removed nodes, keeping the rest, and invalidating the expression's cached hash.

**List-specific.** The rewrite fires only for the head `List`; for any other
head `Nothing` is an ordinary symbol. This is what makes
`Table[If[test, val, Nothing], ...]` the idiom for conditionally building a list:
the iterator produces a plain `List` whose unwanted slots hold `Nothing`, and the
evaluator strips them on the next step.

**Data structures / limits.** O(n) over the list's arguments, one reallocation
when a `Nothing` is present (none otherwise). `Nothing` carries `ATTR_PROTECTED`
(`src/attr.c`). No numeric surface — it is a structural identity element.
