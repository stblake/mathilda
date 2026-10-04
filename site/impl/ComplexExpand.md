---
source: src/complex_expand.c
---
**Algorithm.** `builtin_complex_expand` first peels any trailing
`TargetFunctions -> ...` option rule off the positionals (choosing the output
basis `{Re, Im}` by default, `{Abs, Arg}`, or `Conjugate`), then builds the
complex-variable list from the optional second positional — symbols or
*patterns*; anything not matching one of them is treated as real. The work is a
single recursive routine `cx_decompose(e, ctx, &re, &im)` that writes two
real-valued expressions with `e == re + I*im`. `Plus` sums the pairs
componentwise; `Times` folds them by complex multiplication (`cx_cmul`); `Power`
goes through `cx_power` — `E^z` and a literal `Exp[z]` use
`Exp[u](Cos[v] + I Sin[v])`, a real base with an integer exponent stays real, a
complex base with an integer exponent is raised by square-and-multiply on the
`(re, im)` pair, a positive-real base with a real exponent stays real, and
everything else uses the polar master formula `r^u Exp[-v*theta]` with
`theta = Arg[base]`. `Log[w]` becomes `(1/2) Log[Expand[u^2+v^2]] + I Arg[w]`;
the twelve circular/hyperbolic heads expand through the standard addition-formula
closed forms in `cx_apply_elem` (Tan/Cot/Tanh/Coth via double-angle rational
forms, Sec/Csc/Sech/Csch as reciprocals), and the twelve inverse heads are
rewritten to logarithmic form in `cx_inverse_to_log` and recursed.
`Re/Im/Abs/Arg/Conjugate/Sign/ReIm` are ordinary nodes; an unknown head whose
arguments are all real is assumed real, otherwise it is returned as
`Re[f[...]] + I Im[f[...]]`. A complex atom decomposes to `Re[z]`/`Im[z]` (or
`Abs[z] Cos[Arg z]` / `Abs[z] Sin[Arg z]` under `{Abs, Arg}`); every other atom
is real. `TargetFunctions -> Conjugate` is a separate path: `cx_conjugate_of`
maps `I -> -I` and each complex variable `z -> Conjugate[z]`, then averages the
expression with its conjugate. The front end threads over `List`, the relational
heads and the logical heads (collapsing `Equal`/`Unequal` to `True`/`False` via
`Simplify` when it can), and each leaf returns `Expand[re + I*im]`.

**Data structures.** Everything is `Expr` trees. Every `cx_*` helper borrows its
arguments and returns a freshly-owned, already-evaluated `Expr` built through the
`mk_*`/`ev` (= `eval_and_free`) shorthand; the builtin never frees `res` (the
evaluator owns it). Complex-variable membership is decided by running the pattern
matcher (`match`) against each supplied pattern, so `x_` and richer patterns work
as the variable spec.

**Complexity / limits.** Dominated by the repeated `eval_and_free`/`Expand` at
every node, so cost grows with the expanded size of the real and imaginary parts
(the `Times` fold is quadratic in the number of complex factors). The head is
purely symbolic and structural: no NDArray or `Compile[]` path, and `Protected`
is its only attribute. `Arg[...]`/`Abs[...]` on a symbolic complex atom are kept
inert, so an expression carrying those is as far as the decomposition reduces.
