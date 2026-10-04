### Worked examples

```mathematica
In[1]:= NotElement[2, Integers]  (* 2 is an integer, so the negation is False *)
```

```mathematica
In[1]:= NotElement[1/2, Integers]  (* one half is not an integer, so True *)
```

```mathematica
In[1]:= NotElement[I, Reals]  (* the imaginary unit is not real *)
```

```mathematica
In[1]:= NotElement[{1/2, 3}, Integers]  (* a list is in Integers only if every element is *)
```

```mathematica
In[1]:= NotElement[x, Reals]  (* an undecided membership stays symbolic *)
```

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
