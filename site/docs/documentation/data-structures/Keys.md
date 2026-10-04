# Keys

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Keys[assoc]`**

Gives a list of the keys of an association (or of a rule or list of rules).

**`Keys[{assoc1, assoc2, ...}]`**

Threads over lists of associations and lists of rules.

**`Keys[assoc, f]`**

Wraps each key: {f\[k1\], f\[k2\], ...}.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Keys[<|"a" -> 1, "b" -> 2|>]
Out[1]= {"a", "b"}

In[2]:= Keys[<|"a" -> 1, "b" -> 2|>, f]
Out[2]= {f["a"], f["b"]}

In[3]:= Keys[{<|"a" -> 1|>, <|"b" -> 2, "c" -> 3|>}]
Out[3]= {{"a"}, {"b", "c"}}
```

### Applications (3)

```mathematica
In[4]:= Keys[<|a -> 1, b -> 2, c -> 3|>]
Out[4]= {a, b, c}
```

Wrap each key as f[k]

```mathematica
In[5]:= Keys[<|a -> 1, b -> 2|>, f]
Out[5]= {f[a], f[b]}
```

A bare list of rules is accepted too

```mathematica
In[6]:= Keys[{a -> 1, b -> 2}]
Out[6]= {a, b}
```

## Implementation notes

**Algorithm.** The base `builtin_keys` (`keys_or_values(res, want_keys=true)`) walks
the entries once and copies each rule's key into a fresh `List`. It also accepts a
single `Rule` (returning its key) or a bare list of rules, for Wolfram parity. The
extended forms live in a wrapper installed by `assoc_ops_init` (`ops_keys` →
`keys_values` → `kv_extract`): `Keys[assoc, f]` wraps each key as `f[k]`, and both
forms thread recursively over nested lists of associations and lists of rules.

**Data structures.** Plain `Expr` trees. One `List` of `expr_copy`'d keys is built;
no hash index is needed because the operation is a straight structural read in
insertion order.

**Complexity / limits.** O(n) in the number of entries. A list element that is not a
rule or association yields the `Keys::invrl` message (or, inside an association,
leaves the call unevaluated).

**Attributes:** `Protected`.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_atomicity.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_atomicity.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

Keys are returned in insertion order. `Keys[assoc, f]` maps `f` over each key, and
both forms thread over a list of associations or a list of rules. The dual builtin
`Values` reads the right-hand sides instead.
