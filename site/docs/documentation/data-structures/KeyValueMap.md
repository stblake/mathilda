# KeyValueMap

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyValueMap[f, assoc]`**

Gives {f\[k1, v1\], f\[k2, v2\], ...}.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= KeyValueMap[f, <|a -> 1, b -> 2|>]
Out[1]= {f[a, 1], f[b, 2]}

In[2]:= KeyValueMap[Plus, <|1 -> 10, 2 -> 20|>]
Out[2]= {11, 22}
```

### Applications (2)

F sees each key and value

```mathematica
In[3]:= KeyValueMap[f, <|a -> 1, b -> 2|>]
Out[3]= {f[a, 1], f[b, 2]}
```

Rebuild rules from key and value

```mathematica
In[4]:= KeyValueMap[#1 -> #2^2 &, <|a -> 2, b -> 3|>]
Out[4]= {a -> 4, b -> 9}
```

## Implementation notes

**Algorithm.** `builtin_keyvaluemap` builds `List[f[k1, v1], f[k2, v2], …]` — one
`f` applied to the key and value of each entry, in order. The applications are
*not* forced inside the builtin; the enclosing evaluator reduces each `f[ki, vi]` as
usual, so an operator like `#1 -> #2^2 &` produces ordinary rules.

**Data structures.** One `List` of freshly constructed `f[k, v]` nodes; no hash index
(it is a straight ordered walk).

**Complexity / limits.** O(n) to construct, plus whatever the evaluator spends
reducing each application. Unlike `KeyMap` (transforms keys, result is an
association) and `Map` over an association (transforms values), `KeyValueMap` sees
both parts and returns a plain `List`.

**Attributes:** `Protected`.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

`f` is applied to each key and value together, and the results form a plain `List`
(not an association). This is the two-argument counterpart of mapping over an
association, which sees only the values.
