# HoldAllComplete

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`HoldAllComplete`**

is an attribute which specifies that all arguments to a function are not to be modified or looked at in any way in the process of evaluation.

<details>
<summary>Notes</summary>

HoldAllComplete prevents argument evaluation, Sequence flattening inside arguments, Unevaluated wrapper stripping, and application of Evaluate. Evaluate cannot override HoldAllComplete.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

```mathematica
In[1]:= SetAttributes[h, HoldAllComplete]
Out[1]= Null

In[2]:= h[1+1]
Out[2]= h[1 + 1]

In[3]:= Attributes[h]
Out[3]= {HoldAllComplete}
```

## Implementation notes

`HoldAllComplete` is an attribute name, not a function. It maps to the bitflag `ATTR_HOLDALLCOMPLETE` (`attr_name_to_flag` / `get_attributes` in `src/attr.c`); when set on a symbol the evaluator holds all arguments and additionally bypasses Sequence flattening, `Unevaluated` stripping, and upvalue lookup for that head.

**Attributes:** none registered.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [HoldComplete](../../expression-information/HoldComplete/), [Unevaluated](../../expression-information/Unevaluated/)

- Source: [`src/attr.c`](https://github.com/stblake/mathilda/blob/main/src/attr.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)

## Notes & additional examples

### Notes

`HoldAllComplete` is the strongest hold attribute: it not only suppresses argument evaluation but also blocks `Sequence` flattening, `Unevaluated` stripping, and `Evaluate`, so the arguments are passed through untouched.
