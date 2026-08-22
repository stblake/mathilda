# Append

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Append[expr, elem] adds elem to the end of expr.`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Append[<|"a" -> 1, "b" -> 2|>, "c" -> 3]
Out[1]= <|"a" -> 1, "b" -> 2, "c" -> 3|>

In[2]:= Append[<|"a" -> 1|>, "a" -> 99]
Out[2]= <|"a" -> 99|>
```

## Implementation notes

`builtin_append` (in `src/core.c`) requires a 2-arg call `Append[expr, elem]` whose first argument is an `EXPR_FUNCTION`. It allocates a fresh argument array one slot longer than `expr`, deep-copies every existing argument plus `elem` into it, and rebuilds a new `EXPR_FUNCTION` with the same head. Works on any head, not just `List`. Returns `NULL` (unevaluated) when the first argument is atomic.

**Attributes:** none registered.

## References

**See also:** [Prepend](../../data-structures/Prepend/), [AssociateTo](../../data-structures/AssociateTo/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
- Tests: [`tests/test_core_algebra.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core_algebra.c)
- Tests: [`tests/test_eval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eval.c)
