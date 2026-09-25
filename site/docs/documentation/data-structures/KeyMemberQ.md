# KeyMemberQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyMemberQ[assoc, patt]`**

Gives True if some key of assoc matches the pattern patt (a literal key is an O(1) probe), else False.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= KeyExistsQ[<|"a" -> 1|>, "a"]
Out[1]= True

In[2]:= KeyFreeQ[<|"a" -> 1|>, "b"]
Out[2]= True

In[3]:= KeyMemberQ[<|1 -> "x", "b" -> "y"|>, _Integer]
Out[3]= True

In[4]:= KeyFreeQ[<|"a" -> 1|>, _Integer]
Out[4]= True

In[5]:= KeyExistsQ[<|"a" -> 1|>, _]
Out[5]= False
```

## Implementation notes

- A pattern-free key keeps the O(1) index probe; a key containing a pattern
  construct (`_`, `x_h`, `Alternatives`, `Except`, `PatternTest`, `Condition`,
  ...) is matched against every key.
- `KeyExistsQ` never treats its key as a pattern: `KeyExistsQ[a, _]` looks for
  the literal key `_`, as in Mathematica.

**Attributes:** `Protected`.

## References

**See also:** [KeyExistsQ](../../data-structures/KeyExistsQ/), [KeyFreeQ](../../data-structures/KeyFreeQ/)

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
