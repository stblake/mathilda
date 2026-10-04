---
source: src/core.c
---
**Algorithm.** `builtin_stringq` (`src/core.c`) is a one-argument type test: it
returns `True` exactly when `res->data.function.args[0]->type == EXPR_STRING`, and
`False` for every other leaf or compound. The empty string `""` is still an
`EXPR_STRING`, so it gives `True`. `StringQ` is not `Listable`, so a list of
strings is tested as a single object (a `List` is not a string) and gives `False`
rather than threading.

**Data structures.** A single tag check on the argument `Expr`; nothing is
allocated on the common path beyond the `True`/`False` symbol returned.

**Complexity / limits.** `O(1)`. Any arity other than one is a malformed call
shape: the builtin routes a `StringQ::argx` diagnostic through
`builtin_arg_error` (so `Quiet[]`/`Check[]` see it) and leaves the call
unevaluated. Attributes `Protected`.
