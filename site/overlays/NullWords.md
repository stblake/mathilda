### Worked examples

```mathematica
In[1]:= Head[NullWords]  (* a bare Read/ReadList option symbol *)
```

```mathematica
In[1]:= Options[ReadList, NullWords]  (* the default is False *)
```

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_nw.txt"]; WriteString[str, "a  b\n"]; Close[str];
In[2]:= ReadList["/tmp/mathilda_nw.txt", Word]  (* default: the empty field is skipped *)
In[3]:= ReadList["/tmp/mathilda_nw.txt", Word, NullWords -> True]  (* now kept as "" *)
```

### Notes

`NullWords` is an option for `Read` and `ReadList` deciding whether an empty word is
assumed between two adjacent word separators. It is an inert option-name symbol — no
builtin, no DownValues — read out of a call's option list by the reader.

The default `False` makes the word reader skip the empty field, so the double space above
reads as `{"a", "b"}`; `NullWords -> True` keeps it, giving `{"a", "", "b"}`. It is one of
the tokenisation options (with `WordSeparators`, `TokenWords`, `RecordSeparators`,
`NullRecords`), all defaulting from `Options[Read]`, so `SetOptions[Read, NullWords ->
True]` would change the global default.
