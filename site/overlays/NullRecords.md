### Worked examples

```mathematica
In[1]:= s = OpenWrite["/tmp/mathilda_null.txt"]; WriteString[s, "a,b,,c"]; Close[s]; ReadList["/tmp/mathilda_null.txt", Record, RecordSeparators -> {","}]  (* default: empty field dropped *)
```

```mathematica
In[1]:= s = OpenWrite["/tmp/mathilda_null.txt"]; WriteString[s, "a,b,,c"]; Close[s]; ReadList["/tmp/mathilda_null.txt", Record, RecordSeparators -> {","}, NullRecords -> True]  (* empty field kept *)
```

### Notes

`NullRecords` is an option for `Read` and `ReadList` specifying whether a null (empty)
record is assumed between two adjacent record separators. The default is `False`.

With the default, `"a,b,,c"` split on `","` collapses the doubled separator and yields
`{"a", "b", "c"}`; with `NullRecords -> True` the empty field between the two commas is
kept, giving `{"a", "b", "", "c"}`. It is an inert, `Protected` keyword with a registered
default (so `SetOptions` works), meaningful only with the `Record` read type and
`RecordSeparators`. The word-level analogue is `NullWords`.
