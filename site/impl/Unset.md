---
source: src/core.c
---
**Algorithm.** `builtin_unset` (`src/core.c`) implements `lhs =.`, removing the
single definition whose left-hand side is `lhs` rather than every rule on a
symbol. It first works out which symbol owns the rule and whether it is an
OwnValue or a DownValue: a bare symbol is an OwnValue on itself; `f[...]` is a
DownValue keyed by the head `f`; and because `f[x_] /; cond =.` parses to
`Unset[Condition[f[x_], cond]]`, a `Condition` wrapper is unwrapped to find the
inner head. A non-assignable left-hand side (e.g. `Unset[5]`) returns `NULL`.

A `Protected` owner is refused with `Unset::wrsym` (mirroring `Set`),
returning `Null`. Otherwise `symtab_remove_matching_rule(name, lhs, own_value)`
deletes exactly the rule whose stored pattern equals `lhs`, leaving the symbol's
other definitions, attributes, and remaining rules intact. The head returns
`Null`.

**Attributes & limits.** `Unset` carries `HoldFirst | Protected`, so the target
is not evaluated to its current value before the rule is located. It takes one
argument.
