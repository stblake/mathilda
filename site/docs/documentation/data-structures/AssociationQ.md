# AssociationQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AssociationQ[expr]`**

Gives True if expr is a valid Association (every entry a Rule or RuleDelayed), else False.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= AssociationQ[<|"a" -> 1|>]
Out[1]= True

In[2]:= AssociationQ[{1, 2, 3}]
Out[2]= False

In[3]:= AssociationQ[Association[f[a -> 1]]]
Out[3]= False
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Rule](../../assignment-and-rules/Rule/), [RuleDelayed](../../assignment-and-rules/RuleDelayed/)

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_atomicity.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_atomicity.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
