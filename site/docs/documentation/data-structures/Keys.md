# Keys

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Keys[assoc]
    Gives a list of the keys of an association (or rules).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Keys[<|"a" -> 1, "b" -> 2|>]
Out[1]= {"a", "b"}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
