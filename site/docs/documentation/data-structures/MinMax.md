# MinMax

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MinMax[list]`**

Gives {Min\[list\], Max\[list\]}. Over an association, uses the values.

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= MinMax[<|"a" -> 3, "b" -> 1, "c" -> 9|>]
Out[1]= {1, 9}
```

## Implementation notes

**Attributes:** `Protected`.

## See also

[Min](../../data-structures/Min/), [Max](../../data-structures/Max/)

## References

- Source: [`src/list/list_init.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_init.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
