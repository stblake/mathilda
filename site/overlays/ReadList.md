### Worked examples

```mathematica
In[1]:= Put[2, 3, 5, 7, 11, "/tmp/mathilda_rl.txt"];  (* Put writes one expression per line *)
In[2]:= ReadList["/tmp/mathilda_rl.txt", Number]  (* read every number into a flat list *)
```

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_rl2.txt"]; WriteString[str, "a 1\nb 2\nc 3\n"]; Close[str];
In[2]:= ReadList["/tmp/mathilda_rl2.txt", {Word, Number}]  (* one of each type per pass -> sublists *)
```

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_rl3.txt"]; WriteString[str, "a,b,c\n1,2,3\n"]; Close[str];
In[2]:= ReadList["/tmp/mathilda_rl3.txt", Word, WordSeparators -> {","}]  (* comma-separated fields *)
```

### Notes

`ReadList["file", type]` reads objects of `type` until end of file into a flat
list; `ReadList["file", {t1, ..., tk}]` reads one of each type per pass and groups
each pass into a `k`-element sublist (handy for columnar data). With no type it
reads every remaining expression. A trailing integer `n` stops after `n` objects or
passes.

Unlike `Get`, which only evaluates each expression for its effect,
`ReadList` **collects** the results, and the non-expression read types let a data
file be parsed field by field. `ReadList` loops the shared reading engine that
`Read` uses, so the read types and the separator options
(`RecordSeparators`, `WordSeparators`, `TokenWords`, `NullRecords`, `NullWords`,
from `Options[ReadList]`) are exactly the same. A named file is opened and closed by
`ReadList`; an already-open `InputStream` is read from its current point and left
open.
