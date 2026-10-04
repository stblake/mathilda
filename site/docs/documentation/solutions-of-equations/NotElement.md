# NotElement

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NotElement[x, dom]`**

The statement that x is not an element of the domain dom -- the negation of Element\[x, dom\].  Decides to True or False when the membership decides, and stays symbolic otherwise.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

2 is an integer, so the negation is False

```mathematica
In[1]:= NotElement[2, Integers]
Out[1]= False
```

One half is not an integer, so True

```mathematica
In[2]:= NotElement[1/2, Integers]
Out[2]= True
```

The imaginary unit is not real

```mathematica
In[3]:= NotElement[I, Reals]
Out[3]= True
```

A list is in Integers only if every element is

```mathematica
In[4]:= NotElement[{1/2, 3}, Integers]
Out[4]= True
```

An undecided membership stays symbolic

```mathematica
In[5]:= NotElement[x, Reals]
Out[5]= NotElement[x, Reals]
```

## Implementation notes

**Algorithm.** `builtin_not_element` is the Boolean negation of `Element`. Given
`NotElement[x, dom]` it builds `Element[x, dom]`, evaluates it, and inverts the
verdict: an `Element` that decides to `True` yields `False`, one that decides to
`False` yields `True`. If `Element` does not resolve to a literal (it returns a
symbolic `Element[...]`, or anything that is not `True`/`False`), the builtin
returns `NULL`, so `NotElement[x, dom]` stays symbolic. All the membership
logic — the recognised domains (`Integers`, `Rationals`, `Reals`, `Complexes`,
`Algebraics`, `Booleans`, ...), and the distribution of a container first argument
(a `List` or `Alternatives`) over its members — lives in `Element`; `NotElement`
only flips the result.

In a logical statement fed to `Reduce`, `LogicalExpand`, or the quantifier engine,
`NotElement` is also understood directly by the DNF translator (`to_dnf` /
`to_dnf_neg` recognise `Element`/`NotElement` as dual literals over a container),
so `!Element[...]` and `NotElement[...]` normalise to the same form.

**Data structures.** Trivial: two `Expr` arguments copied into an `Element` node
that is handed to `evaluate`; the result is a single symbol or `NULL`. No
intermediate representation of its own.

**Complexity / limits.** As cheap as the underlying `Element` call. It decides
precisely when `Element` decides (e.g. `NotElement[I, Reals] -> True`,
`NotElement[2, Integers] -> False`) and is left unevaluated otherwise
(`NotElement[x, Reals]`); it performs no reasoning `Element` itself does not.

**Attributes:** `Protected`.

## References

**See also:** [LogicalExpand](../../solutions-of-equations/LogicalExpand/)

- Source: [`src/solve/reduce_companions.c`](https://github.com/stblake/mathilda/blob/main/src/solve/reduce_companions.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)
- Tests: [`tests/test_deriv.c`](https://github.com/stblake/mathilda/blob/main/tests/test_deriv.c)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)

## Notes & additional examples

### Notes

`NotElement[x, dom]` is the statement that `x` is **not** an element of the domain
`dom` — the negation of `Element[x, dom]`. It decides to `True` or `False`
whenever the membership itself decides (`NotElement[I, Reals] -> True`,
`NotElement[3, Reals] -> False`) and otherwise stays symbolic
(`NotElement[x, Reals]`).

All the membership logic lives in `Element`: the recognised domains (`Integers`,
`Rationals`, `Reals`, `Complexes`, `Algebraics`, `Booleans`, ...) and the
distribution of a container first argument — a `List` such as `{1/2, 3}`, or an
`Alternatives` — over its members. `NotElement` simply inverts that verdict.

`NotElement` is also the head `LogicalExpand` emits for a negated membership, so
`!Element[x, dom]` and `NotElement[x, dom]` normalise to the same literal inside
`Reduce`, `LogicalExpand`, and the quantifier engine.
