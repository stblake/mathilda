# HoldFirst

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`HoldFirst`**

is an attribute that specifies that the first argument to a function is to be maintained in an unevaluated form.

<details>
<summary>Notes</summary>

Use Evaluate to evaluate a held first argument in a controlled way.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

An attribute symbol carries itself, plus Protected

```mathematica
In[1]:= Attributes[HoldFirst]
Out[1]= {HoldFirst, Protected}
```

Give a test function the attribute

```mathematica
In[2]:= SetAttributes[holdfn, HoldFirst]
```

The first argument is held, the rest evaluate

```mathematica
In[3]:= holdfn[1 + 1, 2 + 2]
Out[3]= holdfn[1 + 1, 4]
```

And the function now reports it

```mathematica
In[4]:= Attributes[holdfn]
Out[4]= {HoldFirst}
```

## Implementation notes

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

**Attributes:** `HoldFirst`, `Protected`.

## References

- Source: [`src/attr.c`](https://github.com/stblake/mathilda/blob/main/src/attr.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`HoldFirst` is an evaluation attribute. Setting it on a symbol `f` makes the
evaluator keep `f`'s first argument unevaluated while evaluating the others, which
is exactly what assignment heads such as `Set`, `AppendTo` and `SetAttributes` rely
on. It is one bit (`ATTR_HOLDFIRST`) in a symbol's attribute mask; `SetAttributes`
/ `ClearAttributes` flip it and `Attributes` reads it back. `Evaluate` can override
a hold locally. It is the first-argument member of the `HoldFirst` / `HoldRest` /
`HoldAll` family.
