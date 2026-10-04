# List

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

{e1, e2, ...} or List\[e1, e2, ...\] represents an ordered list of the elements ei. List is the fundamental container head: vectors are lists, matrices are lists of lists, and the structural operators (Part, Map, Take, Drop, Length, ...) act on List. Elements are evaluated normally and kept in the given order (List has no Orderless attribute). The parser writes the {...} syntax to List, and the printer renders List\[...\] back as {...}.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

{...} is syntax for List[...]

```mathematica
In[1]:= FullForm[{p, q, r}]
Out[1]= List[p, q, r]
```

```mathematica
In[2]:= Head[{1, 2, 3}]
Out[2]= List
```

The long form and the brace form are identical

```mathematica
In[3]:= List[1, 2, 3] === {1, 2, 3}
Out[3]= True
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Dot 6x6 x 6x6 x 10000 | 338 s | 6.36 s | 4.01 s |
| Sort 4x10^6 | 42.2 s | 68.7 s | 111 s |
| gather v[[idx]], 4x10^6 | 16.8 s | 6.66 s | 7.18 s |
| Union of 4x10^6 integers | 12.4 s | 71.1 s | 376 s |
| Reverse 4x10^6 | 5.37 s | 0.297 s | 0.982 s |
| Inverse 3x3 x 5000 | 2.7 s | 4.75 s | 7.67 s |

## Implementation notes

**Definition.** `List[e1, e2, ...]`, written `{e1, e2, ...}`, is the fundamental
ordered-container head. It is pure structure: there is no `List` builtin — the evaluator
simply evaluates each element and keeps it. Vectors are lists, matrices are lists of
lists, and the structural operators (`Part`, `Map`, `Take`, `Drop`, `Length`, ...) are
defined to act on `List`.

**Representation.** A bare `EXPR_SYMBOL` head (interned `SYM_List`) on an
`EXPR_FUNCTION` node. The parser (`src/parse.c`) rewrites the `{...}` surface syntax to a
`List[...]` node, and the printer (`src/print.c`) renders `List[...]` back as `{...}`;
`FullForm[{a, b, c}]` exposes the underlying `List[a, b, c]`. The packed-array
substrate (`src/pack.c`) can store a numeric `List` as a dense machine buffer
transparently, but the logical head is always `List`.

**Usage & limits.** Elements are evaluated normally and kept in the given order — `List`
has **no** `Orderless` attribute, so `{3, 1, 2}` is not sorted. In fact `Attributes[List]`
is empty (`{}`): it carries no attributes at all. Its uniformity is the point —
everything is an expression, so one set of generic tools works on every list.

**Attributes:** none registered.

## References

- Source: [`src/parse.c`](https://github.com/stblake/mathilda/blob/main/src/parse.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_array_flatten.c`](https://github.com/stblake/mathilda/blob/main/tests/test_array_flatten.c)
- Tests: [`tests/test_backtrack.c`](https://github.com/stblake/mathilda/blob/main/tests/test_backtrack.c)
- Tests: [`tests/test_blas.c`](https://github.com/stblake/mathilda/blob/main/tests/test_blas.c)
- Tests: [`tests/test_coefficient_rules.c`](https://github.com/stblake/mathilda/blob/main/tests/test_coefficient_rules.c)

## Notes & additional examples

### Notes

`List[e1, e2, ...]`, written `{e1, e2, ...}`, is the fundamental ordered container.
Vectors are lists, matrices are lists of lists, and the structural operators (`Part`,
`Map`, `Take`, `Drop`, `Length`, ...) all act on `List`. The parser rewrites `{...}` to
`List[...]` and the printer renders it back, so `FullForm` is where the `List` head
shows through.

Elements are evaluated normally and kept in the given order — `List` has **no**
`Orderless` attribute, so `{3, 1, 2}` stays `{3, 1, 2}`. In fact it carries no
attributes at all (`Attributes[List]` is `{}`). This uniformity is the design: because
everything is an expression, one set of generic tools works on every list.
