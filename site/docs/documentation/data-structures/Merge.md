# Merge

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Merge[{assoc1, assoc2, ...}, f]
    Combines associations, applying f to the list of values collected
    for each key (e.g. Merge[{...}, Total]).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Merge[{<|"a" -> 1|>, <|"a" -> 2, "b" -> 3|>}, Total]
Out[1]= <|"a" -> 3, "b" -> 3|>

In[2]:= Merge[<|"g1" -> <|"a" -> 1|>, "g2" -> <|"a" -> 2, "b" -> 3|>|>, Total]
Out[2]= <|"a" -> 3, "b" -> 3|>
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
