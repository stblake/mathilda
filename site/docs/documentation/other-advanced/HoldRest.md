# HoldRest

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`HoldRest`**

is an attribute that specifies that all but the first argument to a function are to be maintained in an unevaluated form.

<details>
<summary>Notes</summary>

Use Evaluate to evaluate a held argument in a controlled way.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

Set the attribute on a head

```mathematica
In[1]:= ClearAll[g]; SetAttributes[g, HoldRest]; Attributes[g]
Out[1]= {HoldRest}
```

First arg evaluates, the rest are held

```mathematica
In[2]:= ClearAll[g]; SetAttributes[g, HoldRest]; g[1 + 1, 2 + 2]
Out[2]= g[2, 2 + 2]
```

If holds its branches

```mathematica
In[3]:= MemberQ[Attributes[If], HoldRest]
Out[3]= True
```

## Implementation notes

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

**Attributes:** `HoldRest`, `Protected`.

## References

- Source: [`src/attr.c`](https://github.com/stblake/mathilda/blob/main/src/attr.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`HoldRest` is an attribute: when a head has it, every argument *except the first* is kept
unevaluated. In `g[1 + 1, 2 + 2]` above, `1 + 1` evaluates to `2` while `2 + 2` is left
as written, giving `g[2, 2 + 2]`.

It is the complement of `HoldFirst` and the partial form of `HoldAll`. Control
constructs such as `If` carry it so a branch does not evaluate until the condition picks
it; `RuleDelayed` carries it too. Use `Evaluate` to force a held argument to evaluate in
a controlled way. Like all attributes it is global to the head and is cleared by
`ClearAll` or `ClearAttributes`.
