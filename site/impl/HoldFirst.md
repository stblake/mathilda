---
source: src/attr.c
---
**Definition.** `HoldFirst` is one of Mathilda's **evaluation attributes**: set on a
symbol `f`, it tells the evaluator to keep the *first* argument of `f[...]`
unevaluated while evaluating the rest. It is a name the attribute machinery
recognises, not a function — it has no builtin and no value of its own, and it is
itself `Protected` (its own attributes are `{HoldFirst, Protected}`).

**Representation.** Internally it is the bit `ATTR_HOLDFIRST` in a symbol's
attribute bitmask. `string_to_attribute("HoldFirst")` (`src/attr.c`) maps the name
to that bit, so `SetAttributes[f, HoldFirst]` and `ClearAttributes` flip it;
`Attributes[f]` reads the mask back as a list of names. The default attribute table
in `attr.c` already carries it for the built-ins that need it (`Set`, `AppendTo`,
`PrependTo`, `MessageName`, `Pattern`, `SetAttributes`, ...). Changing it bumps the
rule epoch so any cached evaluation of `f` is invalidated.

**Usage & limits.** When the evaluator processes `f[a1, a2, ...]` it reads `f`'s
attributes before evaluating arguments; with `HoldFirst`, `a1` is left as written
and `a2, ...` evaluate normally (step 3 of the evaluator loop). `Evaluate[a1]`
overrides the hold locally. It is the single-argument-held member of the
`HoldFirst` / `HoldRest` / `HoldAll` family; use `HoldAll` to hold every argument,
`HoldRest` to hold all but the first.
