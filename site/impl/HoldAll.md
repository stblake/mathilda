---
source: src/attr.c
---
**What it is.** `HoldAll` is an attribute *symbol*, not a function — one of the tokens the
attribute system (`src/attr.c`) maps to and from bitflags. In `attr.h` the flag is
`ATTR_HOLDALL`, defined as the pair `ATTR_HOLDFIRST | ATTR_HOLDREST`.
`get_attribute_flag("HoldAll")` returns that combined flag, and when
`attributes_to_list` finds both bits set it emits the single symbol `HoldAll` (otherwise
`HoldFirst` / `HoldRest` individually).

**Effect.** When a head carries `ATTR_HOLDALL` the evaluator suppresses evaluation of
every argument before calling the head (step 3 of the evaluation loop), so the arguments
reach the head unevaluated. `Hold`, `SetDelayed`, `Function`, `Attributes`,
`OwnValues`/`DownValues` and many control-flow heads carry it. Unlike `HoldAllComplete`
(`ATTR_HOLDALLCOMPLETE`) it does *not* also suppress upvalue lookup or the
`Sequence`/`Unevaluated`-stripping machinery.

**Usage.** `HoldAll` has no C handler and no standalone meaning; it appears only inside
`Attributes[...]`, `SetAttributes[sym, HoldAll]` and `ClearAttributes[sym, HoldAll]`.
Setting it on a head `f` is what makes `f[2 + 3]` stay `f[2 + 3]` rather than reducing to
`f[5]`.
