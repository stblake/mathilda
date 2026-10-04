# HoldAll

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HoldAll`**

is an attribute that specifies that all arguments to a function are to be maintained in an unevaluated form.

<details>
<summary>Notes</summary>

You can use Evaluate to evaluate the arguments of a HoldAll function in a controlled way. Even when a function has attribute HoldAll, Sequence objects that appear in its arguments are still by default flattened; Unevaluated wrappers on a held argument are, however, left intact (use HoldComplete/HoldAllComplete to also suppress Sequence flattening).

</details>

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Hold[1+1, 2+2]
Out[1]= Hold[1 + 1, 2 + 2]

In[2]:= Hold[Sequence[a, b], c]
Out[2]= Hold[a, b, c]

In[3]:= Hold[Unevaluated[1+2]]
Out[3]= Hold[Unevaluated[1 + 2]]

In[4]:= Hold[Evaluate[1+2], 3+4]
Out[4]= Hold[3, 3 + 4]
```

### Applications (5)

Give f the HoldAll attribute

```mathematica
In[5]:= SetAttributes[f, HoldAll]
```

```mathematica
In[6]:= Attributes[f]
Out[6]= {HoldAll}
```

The argument is held, so the sum is not computed

```mathematica
In[7]:= f[2 + 3]
Out[7]= f[2 + 3]
```

An ordinary head evaluates its argument first

```mathematica
In[8]:= g[2 + 3]
Out[8]= g[5]
```

:= is defined to hold all of its arguments

```mathematica
In[9]:= Attributes[SetDelayed]
Out[9]= {HoldAll, Protected, SequenceHold}
```

## Implementation notes

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

**Attributes:** none registered.

## References

**See also:** [SequenceHold](../../expression-information/SequenceHold/), [HoldFirst](../../other-advanced/HoldFirst/), [HoldRest](../../other-advanced/HoldRest/), [Hold](../../expression-information/Hold/), [HoldForm](../../expression-information/HoldForm/), [HoldPattern](../../pattern-matching/HoldPattern/), [Function](../../functional-programming/Function/)

- Source: [`src/attr.c`](https://github.com/stblake/mathilda/blob/main/src/attr.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)

## Notes & additional examples

### Notes

`HoldAll` is an attribute symbol, not a function. A head carrying it has all of its
arguments held unevaluated before the head runs — the mechanism behind `Hold`,
`SetDelayed` (`:=`), `Function` and the control-flow heads. It appears only inside
`Attributes[...]`, `SetAttributes[sym, HoldAll]` and `ClearAttributes[sym, HoldAll]`.

Internally `HoldAll` is the pair `HoldFirst | HoldRest`; `Attributes` reports that pair as
the single token `HoldAll`. Unlike `HoldAllComplete` it still allows upvalue lookup and
`Sequence`/`Unevaluated` processing on the held arguments.
