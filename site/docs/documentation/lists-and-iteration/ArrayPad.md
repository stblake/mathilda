# ArrayPad

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ArrayPad[array, m]`**

pads array with m elements of 0 on every side of every level.

**`ArrayPad[array, {m, n}]`**

pads with m elements at the start and n at the end of each dimension.

**`ArrayPad[array, {{m1, n1}, {m2, n2}, ...}]`**

pads with mi, ni elements at level i; a negative amount removes elements.

**`ArrayPad[array, amounts, padding]`**

uses the given padding: a constant c, a cyclic list {c1, c2, ...}, or one of "Fixed", "Periodic", "Reflected", "Reversed", "ReversedNegation", "ReflectedDifferences", "ReversedDifferences", "Extrapolated" (which takes the option InterpolationOrder).

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= ArrayPad[{1, 2, 3}, 1]
Out[1]= {0, 1, 2, 3, 0}

In[2]:= ArrayPad[{1, 2, 3}, 2, "Fixed"]
Out[2]= {1, 1, 1, 2, 3, 3, 3}

In[3]:= ArrayPad[{a, b, c}, 3, "Extrapolated"]
Out[3]= {4 a - 3 b, 3 a - 2 b, 2 a - b, a, b, c, -b + 2 c, -2 b + 3 c, -3 b + 4 c}

In[4]:= ArrayPad[Range[10], -2]
Out[4]= {3, 4, 5, 6, 7, 8}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [NDArray](../../linear-algebra/NDArray/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_array_pad.c`](https://github.com/stblake/mathilda/blob/main/tests/test_array_pad.c)
