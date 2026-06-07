---
source: src/core.c
---
**Algorithm.** `builtin_compoundexpression` (in `src/core.c`) walks its argument list in order, calling `evaluate()` on each `expr;` clause and freeing the previous result before storing the next, so only the value of the last clause survives. An empty `CompoundExpression[]` returns the symbol `Null`. The head carries `ATTR_HOLDALL` (registered in `src/attr.c`), so the evaluator does not pre-evaluate the arguments — the builtin evaluates them itself, one at a time, which is what makes `a; b; c` sequence side effects (e.g. assignments installed by `Set`) in source order rather than all at once. After each evaluation the result head is inspected for the control-flow markers `Return`, `Break`, `Continue`, `Throw`, `Abort`, and `Quit`; if any is seen the loop breaks early and returns that wrapper unevaluated so an enclosing construct can act on it.

**Data structures.** Plain iteration over `res->data.function.args`; a single `last_val` accumulator holds the running result.
