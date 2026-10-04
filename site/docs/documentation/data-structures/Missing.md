# Missing

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Missing[]  |  Missing["reason"]  |  Missing["reason", data]`**

Represents missing data. Lookup and key access give Missing\["KeyAbsent", k\] for an absent key; KeyUnion and JoinAcross fill gaps with Missing\[...\]. Test with MissingQ, remove with DeleteMissing.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= <|"a" -> 1|>["q"]
Out[1]= Missing["KeyAbsent", "q"]

In[2]:= Missing["NotAvailable"]
Out[2]= Missing["NotAvailable"]

In[3]:= DeleteMissing[{1, Missing["NotAvailable"], 3}]
Out[3]= {1, 3}
```

### Applications (3)

An inert marker: stays as written

```mathematica
In[4]:= Missing["KeyAbsent", x]
Out[4]= Missing["KeyAbsent", x]
```

An absent key produces one

```mathematica
In[5]:= Lookup[<|a -> 1|>, b]
Out[5]= Missing["KeyAbsent", b]
```

Test with MissingQ

```mathematica
In[6]:= MissingQ[Missing["Unknown"]]
Out[6]= True
```

## Implementation notes

**Algorithm.** `Missing` is an inert, `Protected` head that *represents* missing data;
it has no evaluation builtin, so `Missing[]`, `Missing["reason"]` and
`Missing["reason", data]` persist unchanged. It is produced by the association
machinery: `Lookup` and key access return `Missing["KeyAbsent", key]` for an absent
key (`make_missing` builds it), `KeyUnion` fills padded gaps with it, and
`JoinAcross` fills unmatched rows with `Missing["Unmatched"]`.

**Data structures.** A plain `Expr` function node with head `Missing`; the second
argument, when present, carries the key or reason that explains the absence.

**Complexity / limits.** None — it is a symbolic marker. Test for it with `MissingQ`
and drop such entries from an association with `DeleteMissing`; `Lookup`'s
third-argument default is the usual way to supply a value in place of one.

**Attributes:** `Protected`.

## References

**See also:** [KeyUnion](../../data-structures/KeyUnion/), [JoinAcross](../../data-structures/JoinAcross/), [MissingQ](../../data-structures/MissingQ/), [DeleteMissing](../../data-structures/DeleteMissing/)

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

`Missing[…]` is a placeholder for data that is not there, not something that
evaluates. Absent-key lookups yield `Missing["KeyAbsent", key]`, and `KeyUnion` /
`JoinAcross` fill gaps with `Missing[…]`. Detect it with `MissingQ`, strip it with
`DeleteMissing`, or avoid it entirely by giving `Lookup` a default.
