# AssociationComap

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AssociationComap[{f1, f2, ...}, x]`**

Gives \<|f1 -\> f1\[x\], f2 -\> f2\[x\], ...|\>.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= AssociationComap[{Min, Max, Length}, {3, 1, 2}]
Out[1]= <|Min -> 1, Max -> 3, Length -> 3|>
```

### Applications (2)

```mathematica
In[2]:= AssociationComap[{f, g}, x]
Out[2]= <|f -> f[x], g -> g[x]|>
```

Several summaries of one dataset, labelled by function

```mathematica
In[3]:= AssociationComap[{Total, Max, Min}, {3, 1, 2}]
Out[3]= <|Total -> 6, Max -> 3, Min -> 1|>
```

## Implementation notes

**Algorithm.** `builtin_associationcomap` is the co-map (reverse of
`AssociationMap`): given a `List` of functions `{f1, f2, ...}` and a value `x`,
it builds `<|f1 -> f1[x], f2 -> f2[x], ...|>`. Each function becomes both a key
and the head of an unevaluated application `fi[x]` (built with `mk_call1`), and
the rule array is canonicalised through `assoc_from_rules`.

**Data structures.** A flat rule array handed to `assoc_from_rules`, whose
transient `KeyIndex` de-duplicates the function keys; `x` is copied once per
function into its application.

**Complexity / limits.** `O(k)` for `k` functions (plus the later cost of
evaluating each `fi[x]`). The first argument must be a `List`, otherwise
`AssociationComap::invl` is emitted through the message funnel and the call is
left unevaluated.

**Attributes:** `Protected`.

## References

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)

## Notes & additional examples

### Notes

`AssociationComap[{f1, f2, ...}, x]` is the reverse of `AssociationMap`: it
applies each function to the *same* value `x` and labels the results by the
functions, giving `<|f1 -> f1[x], f2 -> f2[x], ...|>`. It is a compact way to
compute several named summaries of one object — for instance `Total`, `Max`, and
`Min` of a list — in a single keyed record. The first argument must be a list of
functions.
