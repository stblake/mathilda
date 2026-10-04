# RecordSeparators

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`RecordSeparators`**

is an option for Read and ReadList giving the list of strings that separate records. The default is {"\r\n", "\n", "\r"}.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

Split records on commas

```mathematica
In[1]:= s = OpenWrite["/tmp/mathilda_rec.txt"]; WriteString[s, "a,b,,c"]; Close[s]; ReadList["/tmp/mathilda_rec.txt", Record, RecordSeparators -> {","}]
Out[1]= {"a", "b", "c"}
```

Keep the empty field

```mathematica
In[2]:= s = OpenWrite["/tmp/mathilda_rec.txt"]; WriteString[s, "a,b,,c"]; Close[s]; ReadList["/tmp/mathilda_rec.txt", Record, RecordSeparators -> {","}, NullRecords -> True]
Out[2]= {"a", "b", "", "c"}
```

## Implementation notes

**Definition.** `RecordSeparators` is an **option** for `Read` and `ReadList` giving the
list of strings that delimit records (for the `Record` read type). It has no builtin and
no value of its own — it is an inert, `Protected` keyword carried in the reader's option
list, with default `{"\r\n", "\n", "\r"}`.

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_RecordSeparators`). Both
`read.c` and `readlist.c` seed their default option set with
`rule_of(SYM_RecordSeparators, {"\r\n", "\n", "\r"})`; a user `RecordSeparators ->
{...}` rule overrides it. The reader's `next_record` scanner splits on any of these
strings. The companion `NullRecords` option decides whether an empty field between two
adjacent separators is kept.

**Usage & limits.** `Protected`; meaningful only as a `Read`/`ReadList` option. Because
it has a registered default, `SetOptions[Read, RecordSeparators -> ...]` works. Related
tokenisation options are `WordSeparators`, `TokenWords`, `NullRecords` and `NullWords`.

**Attributes:** none registered.

## References

- Source: [`src/io/read.c`](https://github.com/stblake/mathilda/blob/main/src/io/read.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`RecordSeparators` is an option for `Read` and `ReadList` giving the strings that
delimit records (used with the `Record` read type). The default is `{"\r\n", "\n",
"\r"}` — one record per line.

It pairs with `NullRecords`: by default an empty field between two adjacent separators is
dropped, so `"a,b,,c"` split on `","` yields `{"a", "b", "c"}`; with `NullRecords ->
True` the empty field is kept, giving `{"a", "b", "", "c"}`. `RecordSeparators` is an
inert, `Protected` keyword with a registered default, so `SetOptions[Read,
RecordSeparators -> ...]` works. The word-level analogues are `WordSeparators` and
`TokenWords`.
