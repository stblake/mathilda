# Append

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Append[expr, elem] adds elem to the end of expr.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Append[<|"a" -> 1, "b" -> 2|>, "c" -> 3]
Out[1]= <|"a" -> 1, "b" -> 2, "c" -> 3|>

In[2]:= Append[<|"a" -> 1, "b" -> 2|>, "a" -> 9]
Out[2]= <|"b" -> 2, "a" -> 9|>

In[3]:= Prepend[<|"a" -> 1, "b" -> 2|>, "b" -> 9]
Out[3]= <|"b" -> 9, "a" -> 1|>
```

## Implementation notes

`builtin_append` (in `src/core.c`) requires a 2-arg call `Append[expr, elem]` whose first argument is an `EXPR_FUNCTION`. It allocates a fresh argument array one slot longer than `expr`, deep-copies every existing argument plus `elem` into it, and rebuilds a new `EXPR_FUNCTION` with the same head. Works on any head, not just `List`. Returns `NULL` (unevaluated) when the first argument is atomic.

- As in Mathematica 15, an existing entry with the same key is **moved**: it is
  removed from its old position and the new value is placed at the end
  (`Append`) or the front (`Prepend`).
- A non-rule element leaves `Append`/`Prepend` unevaluated; `AppendTo[a, 5]`
  assigns the unevaluated `Append[a, 5]`, exactly like Mathematica.
- `O(n + m)`: membership is tested through an indexed association of the new
  rules.

**Attributes:** none registered.

## References

**See also:** [Prepend](../../data-structures/Prepend/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
