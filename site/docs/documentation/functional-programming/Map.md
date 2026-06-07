# Map

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
f /@ expr or Map[f, expr]
    applies f to each element at level 1 of expr, preserving expr's head.
Map[f, expr, levelspec]
    applies f at the parts of expr selected by levelspec (e.g. {2} for
    level 2 only, Infinity for every level).
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Harold Abelson and Gerald Jay Sussman, *Structure and Interpretation of Computer Programs*, 2nd ed., §2.2.1 (sequence mapping).
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Map[f, {a, b, c}]
Out[1]= {f[a], f[b], f[c]}
```

```mathematica
In[1]:= #^2 & /@ {1, 2, 3, 4}
Out[1]= {1, 4, 9, 16}
```

```mathematica
In[1]:= Map[Reverse, {{1, 2}, {3, 4}}]
Out[1]= {{2, 1}, {4, 3}}
```

```mathematica
In[1]:= Map[f, {{a}, {b}}, {2}]
Out[1]= {{f[a]}, {f[b]}}
```

### Notes

`f /@ expr` is the operator shorthand for `Map[f, expr]`. By default the function
is applied at level 1, i.e. to the immediate elements; a level specification such
as `{2}` reaches deeper into nested lists. `Map` works on any expression, not
only `List` — the head is preserved while each argument is wrapped by `f`. Pure
functions (`#^2 &`) are the idiomatic first argument.
