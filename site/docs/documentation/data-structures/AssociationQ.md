# AssociationQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AssociationQ[expr]`**

Gives True if expr is a valid Association (every entry a Rule or RuleDelayed), else False.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= AssociationQ[<|"a" -> 1|>]
Out[1]= True

In[2]:= AssociationQ[{1, 2, 3}]
Out[2]= False

In[3]:= AssociationQ[Association[f[a -> 1]]]
Out[3]= False
```

### Applications (2)

```mathematica
In[4]:= AssociationQ[<|a -> 1, b -> 2|>]
Out[4]= True
```

A list of rules is not an association

```mathematica
In[5]:= AssociationQ[{a -> 1, b -> 2}]
Out[5]= False
```

## Implementation notes

**Algorithm.** `builtin_associationq` takes a single argument and returns the
symbol `True` or `False` according to `assoc_is_wellformed` (`assoc_struct.h`).
Having the `Association` head is not sufficient: the constructor leaves a
malformed `Association[1, 2]`, or a `Map`/`Apply` result such as
`Association[f[a -> 1]]`, unevaluated, and such a node is an ordinary expression,
not an association. The well-formedness test therefore requires every entry to be
a two-argument `Rule` or `RuleDelayed`.

**Data structures.** A pure structural scan over the argument's child array; no
hash index is built and the argument is only read, never copied.

**Complexity / limits.** `O(n)` in the number of entries (a bare `<|...|>`
literal that already carries a validated shape still re-scans). A bare list of
rules `{a -> 1}` answers `False` — it is a `List`, not an `Association`. The
head is a `*Q` predicate: it always returns a Boolean and never stays
unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [Rule](../../assignment-and-rules/Rule/), [RuleDelayed](../../assignment-and-rules/RuleDelayed/)

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_atomicity.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_atomicity.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)

## Notes & additional examples

### Notes

`AssociationQ[expr]` returns `True` only for a well-formed association — one with
the `Association` head whose every entry is a two-argument `Rule` or
`RuleDelayed`. A bare list of rules such as `{a -> 1, b -> 2}` is a `List`, not an
association, so it answers `False`; likewise a malformed `Association[1, 2]` left
unevaluated by the constructor. Like every `*Q` predicate it always returns a
Boolean, never staying unevaluated.
