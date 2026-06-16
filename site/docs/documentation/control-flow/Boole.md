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

`builtin_boole` is a one-argument map from boolean symbols to integers: it returns `1` for the interned `True`, `0` for `False`, and `NULL` (stays unevaluated) for anything else. It is registered `ATTR_LISTABLE | ATTR_PROTECTED`, so the evaluator threads it over `List` arguments automatically and the handler only needs the scalar case.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/boolean.c`](https://github.com/stblake/mathilda/blob/main/src/boolean.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Boole[3 > 2]
Out[1]= 1
```

```mathematica
In[1]:= Boole[{True, False, True}]
Out[1]= {1, 0, 1}
```

```mathematica
In[1]:= Table[Boole[PrimeQ[n]], {n, 1, 12}]
Out[1]= {0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0}
```

```mathematica
In[1]:= Sum[Boole[GCD[k, 10] == 1], {k, 1, 10}]
Out[1]= 4
```

```mathematica
In[1]:= Sum[Boole[Mod[k, 3] == 0] k^2, {k, 1, 10}]
Out[1]= 126
```

### Notes

`Boole[expr]` is the Iverson bracket: it yields `1` when `expr` is `True` and
`0` when it is `False`, and stays unevaluated otherwise. Being `Listable`, it
threads element-wise over lists. Combined with `Sum` it turns logical predicates
into counting and conditional-summation devices: `Sum[Boole[GCD[k, 10] == 1], {k, 1, 10}]`
counts the integers up to 10 coprime to 10 (Euler's totient phi(10) = 4), and
weighting the bracket by `k^2` sums squares restricted to a residue class.
