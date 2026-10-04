# DeleteMissing

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DeleteMissing[expr]`**

Removes all Missing\[...\] elements (equivalent to DeleteCases\[expr, \_Missing\]). Over an association, drops entries whose value is Missing\[...\].

**`DeleteMissing[expr, n]`**

Removes Missing\[...\] elements at levels 1 through n (n may be Infinity); association values count as one level down.

**`DeleteMissing[expr, n, d]`**

Removes the elements at levels 1..n that contain a Missing\[...\] at depth d or less.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= DeleteMissing[Lookup[<|"a" -> 1, "b" -> 2|>, {"a", "z", "b"}]]
Out[1]= {1, 2}

In[2]:= DeleteMissing[{1, {Missing[], 2}}, 2]
Out[2]= {1, {2}}

In[3]:= DeleteMissing[{{1, Missing[]}, {2}}, 1, 1]
Out[3]= {{2}}
```

### Applications (3)

```mathematica
In[4]:= DeleteMissing[{1, Missing[], 2, Missing["x"], 3}]
Out[4]= {1, 2, 3}

In[5]:= DeleteMissing[{a, Missing[], b, c, Missing["NotAvailable"]}]
Out[5]= {a, b, c}
```

Drops entries whose value is Missing

```mathematica
In[6]:= DeleteMissing[<|a -> 1, b -> Missing[], c -> 3|>]
Out[6]= <|a -> 1, c -> 3|>
```

## Implementation notes

**Algorithm.** `DeleteMissing` drops `Missing[...]` elements. The one-argument
form `DeleteMissing[expr]` (the original `builtin_delete_missing`, in
`patterns.c`) rewrites to `DeleteCases[expr, _Missing]`. The two- and
three-argument forms are handled by `ops_deletemissing` (`assoc_ops.c`):
`DeleteMissing[expr, n]` recurses with `dm_rec` to levels `1..n` (association
values counting one level down), and `DeleteMissing[expr, n, d]` deletes the
level-`1..n` elements that *contain* a `Missing[...]` at depth `d` or less.

**Data structures.** `dm_rec` rebuilds lists and associations bottom-up, dropping
a child when `has_missing_within` finds a qualifying `Missing`; an association
child is re-wrapped with `assoc_entry_with_value` so keys survive. A packed-list
argument is returned unchanged — a machine buffer holds no `Missing`.

**Complexity / limits.** Linear in the tree size visited (bounded by the level
spec `n`). An argument that is not a list, association, or packed list raises
`DeleteMissing::invrp`; a bad level or depth raises `DeleteMissing::arg2` — both
through the message funnel.

**Attributes:** `Protected`.

## References

**See also:** [Lookup](../../data-structures/Lookup/)

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

`DeleteMissing[expr]` removes every `Missing[...]` element — equivalent to
`DeleteCases[expr, _Missing]` — which is the standard way to clean a dataset of
absent-value markers. Over an association it drops the entries whose *value* is
missing, keeping the rest. The two- and three-argument forms
`DeleteMissing[expr, n]` and `DeleteMissing[expr, n, d]` restrict the work to
levels `1..n` and, with `d`, to elements that contain a missing value no deeper
than `d`.
