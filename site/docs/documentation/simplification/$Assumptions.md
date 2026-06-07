# $Assumptions

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
$Assumptions
    is the default setting for the Assumptions option used in Simplify and other functions that take assumptions.
$Assumptions defaults to True (no assumptions). Functions like Assuming temporarily extend $Assumptions for the duration of their body.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= $Assumptions
Out[1]= True
```

## Implementation notes

- A system symbol with default `OwnValue` `True` (no assumptions). `Assuming`

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)
