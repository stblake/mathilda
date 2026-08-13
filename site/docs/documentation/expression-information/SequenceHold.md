# SequenceHold

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`SequenceHold`**

is an attribute which specifies that Sequence objects appearing in the arguments of a function should not automatically be flattened out. The attribute HoldAllComplete implies SequenceHold.

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Sequence](../../expression-information/Sequence/), [HoldAll](../../expression-information/HoldAll/), [HoldAllComplete](../../expression-information/HoldAllComplete/), [Set](../../assignment-and-rules/Set/), [SetDelayed](../../assignment-and-rules/SetDelayed/), [Rule](../../assignment-and-rules/Rule/), [RuleDelayed](../../assignment-and-rules/RuleDelayed/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
