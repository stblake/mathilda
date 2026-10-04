# SequenceHold

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SequenceHold`**

is an attribute which specifies that Sequence objects appearing in the arguments of a function should not automatically be flattened out. The attribute HoldAllComplete implies SequenceHold.

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Hold[Sequence[a, b]]
Out[1]= Hold[a, b]

In[2]:= SetAttributes[h, {HoldAll, SequenceHold}]; h[a, Sequence[b, c]]
Out[2]= h[a, Sequence[b, c]]

In[3]:= splice[x_] := Sequence[x, x, x]; {a, splice[b], c}
Out[3]= {a, b, b, b, c}
```

### Applications (4)

HoldAll alone does not stop the splice

```mathematica
In[4]:= Hold[Sequence[a, b]]
Out[4]= Hold[a, b]
```

Now it is held

```mathematica
In[5]:= SetAttributes[h, {HoldAll, SequenceHold}]; h[a, Sequence[b, c]]
Out[5]= h[a, Sequence[b, c]]
```

Assignment heads carry it

```mathematica
In[6]:= MemberQ[Attributes[Set], SequenceHold]
Out[6]= True
```

Splices at the call site

```mathematica
In[7]:= splice[x_] := Sequence[x, x, x]; {a, splice[b], c}
Out[7]= {a, b, b, b, c}
```

## Implementation notes

**An attribute symbol, not a function.** `SequenceHold` is one of the evaluation
attributes, carried in the `SymbolDef` bitflags as `ATTR_SEQUENCEHOLD`
(`1 << 13`, `src/attr.h`). It names itself in the static `builtin_attrs[]` table
(`src/attr.c`) as a `Protected` symbol, and `attr_from_symbol` maps the symbol
`SequenceHold` back to the bit, so `SetAttributes[h, SequenceHold]` and
`Attributes[h]` round-trip through it.

**Effect.** The evaluator's `Sequence`-splicing pass (step 2.5, `src/eval.c`) is
guarded by `if (hold_all_complete || (attrs & ATTR_SEQUENCEHOLD))`: when the head
being evaluated carries the bit, `Sequence[...]` objects in its arguments are left
unflattened. `HoldAll` alone does **not** suppress splicing — `Hold[Sequence[a,
b]]` still gives `Hold[a, b]` — so `SequenceHold` is the specific lever that stops
it. `HoldAllComplete` implies it (the same gate). The assignment and rule heads
`Set`, `SetDelayed`, `Rule`, and `RuleDelayed` carry `SequenceHold` in the table,
which is what lets a definition store and later deliver a spliceable `Sequence`.

**Scope.** The bit only changes the splicing decision; it is independent of the
`Hold*` family, which governs argument *evaluation*. Any user head that sets it is
honoured automatically by the same gate — there is no per-head special-casing.
Attributes of the symbol itself: `Protected`.

**Attributes:** `Protected`.

## References

**See also:** [Sequence](../../expression-information/Sequence/), [HoldAll](../../expression-information/HoldAll/), [HoldAllComplete](../../expression-information/HoldAllComplete/), [Set](../../assignment-and-rules/Set/), [SetDelayed](../../assignment-and-rules/SetDelayed/), [Rule](../../assignment-and-rules/Rule/), [RuleDelayed](../../assignment-and-rules/RuleDelayed/)

- Source: [`src/attr.c`](https://github.com/stblake/mathilda/blob/main/src/attr.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)

## Notes & additional examples

### Notes

`SequenceHold` is an **attribute**, not a function: it tells the evaluator not to
flatten `Sequence[...]` objects that appear among a function's arguments.
`HoldAll` on its own does *not* prevent the splice — `Hold[Sequence[a, b]]` gives
`Hold[a, b]` — so `SequenceHold` is the specific lever that stops it, and
`HoldAllComplete` implies it.

The assignment and rule heads `Set`, `SetDelayed`, `Rule`, and `RuleDelayed` all
carry `SequenceHold`. That is what lets a definition store a `Sequence` on its
right-hand side and have it splice only when the defined symbol is later used — the
`splice[x_] := Sequence[x, x, x]` idiom. Any user head that sets the attribute is
honoured the same way.
