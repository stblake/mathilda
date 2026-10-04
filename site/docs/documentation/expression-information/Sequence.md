# Sequence

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Sequence[e1, e2, ...]`**

represents a sequence of arguments that is automatically spliced into the argument list of any enclosing function. Sequence\[\] evaporates and Sequence\[e\] acts like the identity. Splicing is suppressed for heads with the attribute SequenceHold or HoldAllComplete.

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= f[a, Sequence[b, c], d]
Out[1]= f[a, b, c, d]

In[2]:= {a, Sequence[b], c, Identity[d]}
Out[2]= {a, b, c, d}

In[3]:= {a, b, g[x, y], h[w], g[z, y]} /. g -> Sequence
Out[3]= {a, b, x, y, h[w], z, y}
```

### Applications (5)

Spliced into the enclosing argument list

```mathematica
In[4]:= f[a, Sequence[b, c], d]
Out[4]= f[a, b, c, d]
```

A one-element Sequence is the identity

```mathematica
In[5]:= {a, Sequence[b], c, Identity[d]}
Out[5]= {a, b, c, d}
```

The empty Sequence evaporates

```mathematica
In[6]:= {a, Sequence[], c}
Out[6]= {a, c}
```

```mathematica
In[7]:= {a, b, g[x, y], h[w], g[z, y]} /. g -> Sequence
Out[7]= {a, b, x, y, h[w], z, y}
```

BlankSequence binds to a Sequence object

```mathematica
In[8]:= f[a, b, c] /. f[x__] -> x
Out[8]= Sequence[a, b, c]
```

## Implementation notes

**Algorithm.** `Sequence` has no C builtin — it is a structural marker the
evaluator splices. `flatten_sequences` (`src/eval.c`) scans a function's argument
list for any child whose head is `SYM_Sequence`; if it finds one it rebuilds the
argument array, copying each `Sequence`'s contents in place of the wrapper, and
frees the emptied wrapper. `Sequence[]` therefore contributes zero arguments (it
evaporates) and `Sequence[e]` contributes one (it acts as the identity). The pass
runs at evaluation step 2.5, **before** `Flat` / `Listable` / `Orderless`, so
`f[a, Sequence[b, c], d]` is `f[a, b, c, d]` by the time attributes are consulted.

**Data structures.** The splice allocates one fresh `Expr**` argument array of the
combined length and re-homes the surviving pointers into it (`expr_copy` for the
spliced children, pointer move for the rest), then invalidates the node's cached
hash. `Sequence` is also the object produced by `BlankSequence` /
`BlankNullSequence` bindings and by `SlotSequence` (`##`), so the same splice
machinery resolves `f[a,b,c] /. f[x__] -> x`.

**Complexity / limits.** `O(n)` in the argument count of the enclosing call.
Splicing is suppressed when the enclosing head carries `SequenceHold` or
`HoldAllComplete` (checked at the step-2.5 gate), which is how an assignment or
rule can carry a `Sequence` that splices only at the eventual call site. A bare
top-level `Sequence[...]` with no enclosing function is left intact. Attributes
`Protected`.

- Attributes: `{Protected}`.
- Splicing happens structurally during evaluation, before `Flat`/`Listable`/
  `Orderless`: `f[a, Sequence[b, c], d]` becomes `f[a, b, c, d]`.
- `Sequence[]` evaporates and `Sequence[e]` acts like the identity, so
  `{a, Sequence[b], c}` gives `{a, b, c}` and `{Sequence[], a}` gives `{a}`.
- A bare `Sequence[...]` with no enclosing function (including one stored in an
  `OwnValue`) is left as a `Sequence` object; it only splices at a call site.
- `Sequence` is the wrapper produced by `BlankSequence`/`BlankNullSequence`
  (`f[a, b, c] /. f[x__] -> x` gives `Sequence[a, b, c]`) and by `SlotSequence`
  (`##& [a, b, c]` gives `Sequence[a, b, c]`).
- Splicing is suppressed for heads carrying `SequenceHold` or `HoldAllComplete`.

**Attributes:** `Protected`.

## References

**See also:** [Flat](../../expression-information/Flat/), [Orderless](../../expression-information/Orderless/), [BlankSequence](../../pattern-matching/BlankSequence/), [BlankNullSequence](../../pattern-matching/BlankNullSequence/), [SlotSequence](../../functional-programming/SlotSequence/), [SequenceHold](../../expression-information/SequenceHold/), [HoldAllComplete](../../expression-information/HoldAllComplete/)

- Source: [`src/eval.c`](https://github.com/stblake/mathilda/blob/main/src/eval.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_assoc_atomicity.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_atomicity.c)
- Tests: [`tests/test_eval_eager_exit.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eval_eager_exit.c)
- Tests: [`tests/test_evaluate.c`](https://github.com/stblake/mathilda/blob/main/tests/test_evaluate.c)
- Tests: [`tests/test_expr_pool.c`](https://github.com/stblake/mathilda/blob/main/tests/test_expr_pool.c)

## Notes & additional examples

### Notes

`Sequence[e1, e2, ...]` is a run of arguments that is automatically spliced into
the argument list of whatever function encloses it. Splicing happens structurally
during evaluation, **before** `Flat` / `Listable` / `Orderless`, so
`f[a, Sequence[b, c], d]` becomes `f[a, b, c, d]`. `Sequence[]` contributes
nothing and `Sequence[e]` contributes one argument.

A bare `Sequence[...]` with no enclosing function — including one stored in an
`OwnValue` — is left as a `Sequence` object and splices only when it reaches a call
site. It is the wrapper produced by `BlankSequence` / `BlankNullSequence` and by
`SlotSequence` (`##`). Splicing is suppressed under a head carrying `SequenceHold`
or `HoldAllComplete`.
