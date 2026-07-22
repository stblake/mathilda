---
source: src/simp/simp_builtins.c
---
**Algorithm.** `builtin_assuming` (`Assuming[assum, body]`, `ATTR_HOLDREST`)
desugars to `Block[{$Assumptions = $Assumptions && assum}, body]` and evaluates
that block. A `List` of assumptions is first normalised to an `And` conjunction
(the standard convention). Building it as `Set[$Assumptions, And[$Assumptions,
assum]]` inside the `Block` variable list reuses Block's existing scope save /
restore machinery, so `$Assumptions` is temporarily extended for the dynamic
extent of `body` and restored afterward. Nested `Assuming` calls compose
naturally because each Block reads the current `$Assumptions` value before
extending it.

**Data structures.** No state of its own — it constructs a `Block[...]` `Expr*`
and hands it to the evaluator; the assumption set lives in the `$Assumptions`
OwnValue.
