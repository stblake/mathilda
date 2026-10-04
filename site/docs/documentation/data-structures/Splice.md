# Splice

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Splice[{e1, e2, ...}]`**

Is replaced by the sequence e1, e2, ... when it appears inside a List or an Association.

**`Splice[{e1, e2, ...}, h]`**

Splices into any head matching the pattern h.

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= {1, Splice[{2, 3}], 4}
Out[1]= {1, 2, 3, 4}

In[2]:= <|"a" -> 1, Splice[{"b" -> 2, "c" -> 3}]|>
Out[2]= <|"a" -> 1, "b" -> 2, "c" -> 3|>

In[3]:= Table[Splice[{i, -i}], {i, 3}]
Out[3]= {1, -1, 2, -2, 3, -3}

In[4]:= f[1, Splice[{2, 3}]]
Out[4]= f[1, Splice[{2, 3}]]

In[5]:= f[1, Splice[{2, 3}, _]]
Out[5]= f[1, 2, 3]
```

### Applications (3)

Splices into the surrounding list

```mathematica
In[6]:= {1, Splice[{2, 3}], 4}
Out[6]= {1, 2, 3, 4}
```

The head pattern _ lets it splice into any head

```mathematica
In[7]:= f[a, Splice[{b, c}, _], d]
Out[7]= f[a, b, c, d]
```

Inert on its own

```mathematica
In[8]:= Splice[{1, 2}]
Out[8]= Splice[{1, 2}]
```

## Implementation notes

**Algorithm.** `Splice` is an evaluator hook, not an ordinary builtin. During the
evaluator's Sequence-flattening step (which runs for heads without `SequenceHold` /
`HoldAllComplete`, and only once a `Splice` argument has actually been seen, so no
other call pays for it), `eval_splice_args` replaces each qualifying `Splice[{e…}]`
with its elements. `splice_applies` gates this: the one-argument `Splice[list]`
splices only inside a `List` or an `Association`; the two-argument `Splice[list, h]`
splices into any head matching the pattern `h` (checked by `MatchQ`), but never into
an association. Anywhere else `Splice` stays inert.

**Data structures.** A `hit` bitmap marks the arguments to expand; the enclosing
node's argument array is rebuilt in place. A packed inner list is materialised to a
nested `List` first.

**Complexity / limits.** O(total arguments after expansion). Because the hook keys on
a `Splice` head actually appearing among the arguments, calls with no `Splice` are
untouched. A bare `Splice[{…}]` at top level (not inside a List/Association, no head
spec) is left unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [Association](../../data-structures/Association/), [SequenceHold](../../expression-information/SequenceHold/), [HoldAllComplete](../../expression-information/HoldAllComplete/)

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)

## Notes & additional examples

### Notes

The one-argument `Splice[{e1, e2, …}]` expands into the enclosing `List` or
`Association`; to splice into another head you give a head pattern, as in
`Splice[list, _]`. On its own, or in a head it does not apply to, `Splice` stays
unevaluated. The work happens in the evaluator's sequence-flattening pass.
