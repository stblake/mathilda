---
source: src/simp/simp.c
---
**Algorithm.** `$Assumptions` is a global assumption store, not a builtin. In
`simp_init` it is given an OwnValue defaulting to `True`
(`symtab_add_own_value("$Assumptions", ...)`). Simplify and Element read it via
`read_dollar_assumptions`, which fetches the OwnValue's replacement *directly*
(`symtab_get_own_values`) and deep-copies it rather than evaluating it — evaluating
would re-fire the rule and, for a bound `Element[x, Reals]`, cause Element to
recurse. `Assuming` extends it by desugaring to `Block[{$Assumptions =
$Assumptions && a}, body]`. The accumulated expression is parsed into an
`AssumeCtx` fact set by `assume_ctx_from_expr` (`simp_assume.c`), which flattens
`And`/`List` conjunctions and splits `Element[{...}, dom]` into per-variable facts.

**Data structures.** A single OwnValue `Rule` on the `$Assumptions` symbol in the
symbol table; consumed as an `AssumeCtx` (flat `Expr*` fact array) during
simplification.
