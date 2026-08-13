# Sequence

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Sequence[e1, e2, ...]`**

represents a sequence of arguments that is automatically spliced into the argument list of any enclosing function. Sequence\[\] evaporates and Sequence\[e\] acts like the identity. Splicing is suppressed for heads with the attribute SequenceHold or HoldAllComplete.

## Examples (3)

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

## Implementation notes

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

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_eval_eager_exit.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eval_eager_exit.c)
- Tests: [`tests/test_evaluate.c`](https://github.com/stblake/mathilda/blob/main/tests/test_evaluate.c)
- Tests: [`tests/test_expr_pool.c`](https://github.com/stblake/mathilda/blob/main/tests/test_expr_pool.c)
- Tests: [`tests/test_expr_sharing.c`](https://github.com/stblake/mathilda/blob/main/tests/test_expr_sharing.c)
