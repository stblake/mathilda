# Values

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Values[assoc]
    Gives a list of the values of an association (or rules).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Values[<|"a" -> 1, "b" -> 2|>]
Out[1]= {1, 2}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
