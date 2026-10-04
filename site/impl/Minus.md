---
source: src/minus.c
---
**Algorithm.** The parser already reads `-x` as `Times[-1, x]`, so `Minus` only
has to exist as a *callable* head. `builtin_minus` returns `Times[-1,
expr_copy(arg)]` and lets `Times` do all the arithmetic — integers, Rationals,
Complex, `Plus` distribution, packed/NDArray buffers. This single-rewrite form
is precisely what makes `Minus` usable as a sort key (`SortBy[list, Minus]`,
`KeySortBy[a, Minus]`). Any argument count other than one emits `Minus::argx`
(through `builtin_arg_error`) and the call is left unevaluated, matching
Wolfram Language, which prints `Minus[x, y]` as `x - y` but does not evaluate it.

**Data structures.** Pure `Expr`: one freshly-built `Times[-1, x]` node, with
the argument deep-copied so the evaluator can still free the original call. Since
the result is simply a `Times`, every surface `Times` already supports is
inherited for free — `Minus` is on `pack.c`'s `AWARE` list for exactly this
reason (a packed or visible `NDArray` argument stays on the buffer once `Times`
sees it), and `Compile[]` lowers the negation at both scalar and rank-1 array
shapes (`compile_emit_arith.c`, with the machine- and GMP-integer domains in
`compile_mgd.c` and the type rule in `compile_infer.c`).

**Complexity / limits.** `O(1)` to build the node; the real cost is whatever
`Times[-1, x]` then does. Attributes are `Listable | NumericFunction |
Protected`, so `Minus` threads over a list and `CompileDiagnostics` reports
`Compiled -> True` at scalar and rank-1 shapes. The same module also defines
`` Internal`SyntacticNegativeQ ``, the leading-minus-sign predicate the ported
integrator's term ordering is built on.
