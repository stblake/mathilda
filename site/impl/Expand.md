---
source: src/expand.c
---
**Algorithm.** `builtin_expand` (in `src/expand.c`) calls `expr_expand_patt(e, patt)`, a structural recursion that distributes products over sums. The dispatch by head:

- **Plus** — expand each summand and rebuild via `eval_and_free`.
- **Times** — expand each factor, then multiply them pairwise with `multiply_all`, a divide-and-conquer fold whose leaf operation `multiply_two` handles the four `Plus×Plus / Plus×atom / atom×Plus / atom×atom` cases, producing every cross term (`a_i · b_j`) and re-summing them.
- **Power[base, k]** with `k` a positive integer `< 100` — expand the base, then `power_expand` raises it by **repeated squaring** (binary exponentiation), each multiply going through the distributing `multiply_two`. (Note: there is no explicit binomial-coefficient formula; the binomial expansion of `(a+b)^k` falls out of the iterated distribution.) Negative or non-integer exponents are left untouched.
- **List / equations / inequalities / And / Or / Not** — threads into each argument (passing operator-symbol slots of `Inequality` through verbatim).

A two-argument `Expand[expr, patt]` only expands subexpressions that *contain* `patt` (gated by `expr_contains_patt`, which uses the pattern matcher); everything else is copied unchanged. Expand does **not** descend into the arguments of arbitrary function heads.

**Data structures.** Plain `Expr*` trees; transient `Expr**` argument buffers per node. All arithmetic rebuilding goes back through the evaluator (`eval_and_free`) so `Plus`/`Times` canonicalisation (Flat/Orderless, like-term collection) applies after each step.

**Complexity / limits.** Worst case is the full cross-product blow-up of distribution: expanding `(x_1+...+x_m)^k` produces `O(m^k)`-sized output, and `Power` is capped at exponent `< 100` to bound it.
