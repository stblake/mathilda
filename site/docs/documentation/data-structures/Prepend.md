# Prepend

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Prepend[expr, elem] adds elem to the beginning of expr.`**

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

**Algorithm.** `builtin_prepend` builds a new function node with the same head whose first
element is a copy of the new element followed by copies of all original arguments. Two-argument
form only; returns `NULL` if the first argument is an atom. (`PrependTo` is the mutating
variant that writes the result back to a symbol's OwnValue.)

- As in Mathematica 15, an existing entry with the same key is **moved**: it is
  removed from its old position and the new value is placed at the end
  (`Append`) or the front (`Prepend`).
- A non-rule element leaves `Append`/`Prepend` unevaluated; `AppendTo[a, 5]`
  assigns the unevaluated `Append[a, 5]`, exactly like Mathematica.
- `O(n + m)`: membership is tested through an indexed association of the new
  rules.

**Attributes:** none registered.

## References

**See also:** [Append](../../data-structures/Append/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_eval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eval.c)
