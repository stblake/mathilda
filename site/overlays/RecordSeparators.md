### Worked examples

```mathematica
In[1]:= s = OpenWrite["/tmp/mathilda_rec.txt"]; WriteString[s, "a,b,,c"]; Close[s]; ReadList["/tmp/mathilda_rec.txt", Record, RecordSeparators -> {","}]  (* split records on commas *)
```

```mathematica
In[1]:= s = OpenWrite["/tmp/mathilda_rec.txt"]; WriteString[s, "a,b,,c"]; Close[s]; ReadList["/tmp/mathilda_rec.txt", Record, RecordSeparators -> {","}, NullRecords -> True]  (* keep the empty field *)
```

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
