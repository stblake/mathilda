# Scan

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Scan[f, expr]`**

Applies f to each element of expr for its side effects and returns Null, discarding the results.

**`Scan[f, expr, levelspec]`**

Applies f to the parts of expr selected by levelspec (default {1}), depth-first with leaves before roots. Option Heads-\>True also visits heads. Throw exits to an enclosing Catch; Return\[ret\] makes the final value ret. Over an association, applies f to each value.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= s = 0; Scan[(s = s + #) &, <|"a" -> 1, "b" -> 2, "c" -> 3|>]; s
Out[1]= 6

In[2]:= Scan[Print, {{{a}}}, Infinity] a {a} {{a}}
Out[2]= {{Null a^3}}

In[3]:= Catch[Scan[If[# > 5, Throw[#]] &, {2, 4, 6, 8}]]
Out[3]= 6
```

### Options (1)

```mathematica
In[4]:= Scan[Print, {a, b}, Heads -> True] List a b
Out[4]= List Null a b
```

## Implementation notes

**Attributes:** `Protected`.

## See also

[Map](../../data-structures/Map/), [Throw](../../control-flow/Throw/), [Catch](../../control-flow/Catch/), [NDArray](../../linear-algebra/NDArray/)

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_catch_throw.c`](https://github.com/stblake/mathilda/blob/main/tests/test_catch_throw.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_core_algebra.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core_algebra.c)
