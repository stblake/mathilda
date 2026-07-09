# KeyMap

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
KeyMap[f, assoc]
    Applies f to every key, keeping values (keys that
    collide collapse with last-value-wins).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= KeyMap[f, <|1 -> 10, 2 -> 20|>]
Out[1]= <|f[1] -> 10, f[2] -> 20|>
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
