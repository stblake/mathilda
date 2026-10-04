# NullRecords

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`NullRecords`**

is an option for Read and ReadList specifying whether a null record is assumed between two adjacent record separators. The default is False.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

Default: empty field dropped

```mathematica
In[1]:= s = OpenWrite["/tmp/mathilda_null.txt"]; WriteString[s, "a,b,,c"]; Close[s]; ReadList["/tmp/mathilda_null.txt", Record, RecordSeparators -> {","}]
Out[1]= {"a", "b", "c"}
```

Empty field kept

```mathematica
In[2]:= s = OpenWrite["/tmp/mathilda_null.txt"]; WriteString[s, "a,b,,c"]; Close[s]; ReadList["/tmp/mathilda_null.txt", Record, RecordSeparators -> {","}, NullRecords -> True]
Out[2]= {"a", "b", "", "c"}
```

## Implementation notes

**Definition.** `NullRecords` is an **option** for `Read` and `ReadList` specifying
whether a null (empty) record is assumed between two adjacent record separators. It has
no builtin and no value of its own — an inert, `Protected` keyword carried in the
reader's option list, with default `False`.

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_NullRecords`). Both `read.c` and
`readlist.c` seed their default option set with `rule_of(SYM_NullRecords, False)`. When
`True`, the `Record` reader keeps the empty field that lies between two consecutive
`RecordSeparators`, so `"a,b,,c"` split on `","` yields `{"a", "b", "", "c"}` rather
than `{"a", "b", "c"}`.

**Usage & limits.** `Protected`; meaningful only as a `Read`/`ReadList` option in
combination with the `Record` read type and `RecordSeparators`. It has a registered
default, so `SetOptions` on it works. The word-level analogue is `NullWords`.

**Attributes:** none registered.

## References

- Source: [`src/io/read.c`](https://github.com/stblake/mathilda/blob/main/src/io/read.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`NullRecords` is an option for `Read` and `ReadList` specifying whether a null (empty)
record is assumed between two adjacent record separators. The default is `False`.

With the default, `"a,b,,c"` split on `","` collapses the doubled separator and yields
`{"a", "b", "c"}`; with `NullRecords -> True` the empty field between the two commas is
kept, giving `{"a", "b", "", "c"}`. It is an inert, `Protected` keyword with a registered
default (so `SetOptions` works), meaningful only with the `Record` read type and
`RecordSeparators`. The word-level analogue is `NullWords`.
