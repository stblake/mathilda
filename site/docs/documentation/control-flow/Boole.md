# Boole

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Boole[expr]
    yields 1 if expr is True and 0 if it is False.
Boole is also known as the Iverson bracket, indicator function, and characteristic function.
Boole is typically used to express integrals and sums over regions given by logical combinations of predicates, and as a dummy-variable encoding for categorical variables in statistics.
Boole[expr] remains unchanged if expr is neither True nor False.
Boole[expr] is effectively equivalent to If[expr, 1, 0].
Boole has attributes {Listable, Protected}, so Boole[{e1, e2, ...}] automatically threads to {Boole[e1], Boole[e2], ...}.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= {Boole[False], Boole[True]}
Out[1]= {0, 1}

In[2]:= Boole[{True, False, True, True, False}]
Out[2]= {1, 0, 1, 1, 0}

In[3]:= Boole[x]
Out[3]= Boole[x]

In[4]:= Total[Boole[# > 0 & /@ {-1, 2, -3, 4, 5}]]
Out[4]= 3
```

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
