---
source: src/attr.c
---
**Definition.** `HoldRest` is an **attribute symbol**. Set on a function head, it tells
the evaluator to keep every argument *except the first* unevaluated. The first argument
is evaluated as usual. It has no builtin — it is a name recognised by the attribute
system and surfaced through `Attributes` / `SetAttributes`.

**Representation.** A bare `EXPR_SYMBOL` that maps to the bitflag `ATTR_HOLDREST` in
`src/attr.c`. During `evaluate_step`, after the head's attributes are read, the
argument-evaluation loop skips the held positions: with `HoldRest` only argument 1 is
evaluated and arguments 2..n are passed through literally. `Hold`-family flags compose
with the others (`Listable`, `Flat`, ...) the generic evaluator consults.

**Usage & limits.** It is the complement of `HoldFirst` and the partial form of
`HoldAll`. Built-ins that carry it include control constructs such as `If` (the branches
must not evaluate until the condition selects one) and `RuleDelayed`. Use `Evaluate` to
force a held argument to evaluate in a controlled way. Like all attributes it is global
to the head and is cleared by `ClearAll` or `ClearAttributes`.
