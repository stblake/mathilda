# Nothing

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Nothing`**

is a symbol that is automatically removed from any list in which it appears as an element: {a, Nothing, b} evaluates to {a, b}. It is the identity element for list construction, so Table\[If\[test, val, Nothing\], ...\] builds a list of just the values for which test held. Any Nothing\[...\] form is removed likewise. Non-list heads treat Nothing as an ordinary symbol.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= {1, Nothing, 2}
Out[1]= {1, 2}

In[2]:= Table[If[EvenQ[i], i, Nothing], {i, 6}]
Out[2]= {2, 4, 6}
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_sequence.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sequence.c)
