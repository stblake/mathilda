---
source: src/funcprog.c
---
**Algorithm.** `builtin_through` implements `Through[p[f1,f2,...][x,y,...]] ->
p[f1[x,y,...], f2[x,y,...], ...]`: it threads the argument tuple of the outer
call through each operator sitting in the (compound) head. For the one-argument
form it requires the head itself to be a function `p[f1,...]`; it builds each
`fi[args...]`, evaluates it, and reassembles under `p`. The optional second
argument `h` restricts the rewrite to heads equal to `h`: `transform_head`
recurses through the head structure, applying the threading only at sub-heads
whose head matches `h` (via `expr_eq`) and copying everything else.

Whenever a transformation actually happened the rebuilt expression is passed
through `evaluate()` so the inner `fi[...]` calls take effect; if nothing matched
the original expression is returned unchanged (`transformed` flag). Atoms are
returned by `expr_copy`.

**Data structures.** `Expr`-tree manipulation only; each inner call is built with
`expr_new_function` over copied argument arrays.
