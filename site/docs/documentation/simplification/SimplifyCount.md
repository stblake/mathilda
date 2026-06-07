# SimplifyCount

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SimplifyCount[expr]
    The complexity measure used by Simplify when no
    ComplexityFunction option (or ComplexityFunction -> Automatic) is
    given. Counts subexpressions; integers contribute their decimal
    digit count plus a constant for the sign. Real numbers contribute
    2 (NumberQ but not Integer/Rational).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= SimplifyCount[100 Log[2]]
Out[1]= 6

In[2]:= SimplifyCount[Log[2^100]]
Out[2]= 32

In[3]:= SimplifyCount[1/2]
Out[3]= 3

In[4]:= SimplifyCount[3.14]
Out[4]= 2
```

## Implementation notes

- `Listable`, `Protected`.
- Per node: a symbol, the integer `0`, or a string counts `1`; a positive integer

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/simp/simp.c`](https://github.com/stblake/mathilda/blob/main/src/simp/simp.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)
