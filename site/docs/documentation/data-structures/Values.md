# Values

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Values[assoc]`**

Gives a list of the values of an association (or of a rule or list of rules).

**`Values[{assoc1, assoc2, ...}]`**

Threads over lists of associations and lists of rules.

**`Values[assoc, f]`**

Wraps each value: {f\[v1\], f\[v2\], ...}.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Values[<|"a" -> 1, "b" -> 2|>]
Out[1]= {1, 2}

In[2]:= Values[<|"a" -> 1, "b" -> 2|>, f]
Out[2]= {f[1], f[2]}

In[3]:= Values[{<|"a" -> 1|>, {"b" -> 2}}]
Out[3]= {{1}, {2}}
```

### Applications (2)

```mathematica
In[4]:= Values[<|a -> 1, b -> 2, c -> 3|>]
Out[4]= {1, 2, 3}
```

Wrap each value as f[v]

```mathematica
In[5]:= Values[<|a -> 1, b -> 2|>, f]
Out[5]= {f[1], f[2]}
```

## Implementation notes

**Algorithm.** The base `builtin_values` (`keys_or_values(res, want_keys=false)`)
walks the entries once and copies each rule's right-hand side into a fresh `List`,
also accepting a single `Rule` or a bare list of rules. The extended forms come from
the `assoc_ops_init` wrapper (`ops_values` → `keys_values` → `kv_extract`):
`Values[assoc, f]` wraps each value as `f[v]`, and both forms thread recursively over
nested lists of associations and lists of rules.

**Data structures.** Plain `Expr` trees; one `List` of `expr_copy`'d values in
insertion order. No hash index — it is a structural read.

**Complexity / limits.** O(n). A non-rule, non-association list element triggers
`Values::invrl` (or leaves the call unevaluated inside an association).

**Attributes:** `Protected`.

## References

**See also:** [Keys](../../data-structures/Keys/)

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)

## Notes & additional examples

### Notes

Values are returned in the association's insertion order, so `Keys` and `Values`
line up positionally. `Values[assoc, f]` maps `f` over each value, and both forms
thread over a list of associations or rules.
