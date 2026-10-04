### Worked examples

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_read.txt"]; WriteString[str, "1 2\n3 4\nx + 1\n"]; Close[str];  (* a small data file *)
In[2]:= ins = OpenRead["/tmp/mathilda_read.txt"];
In[3]:= Read[ins, {{Number, Number}, {Number, Number}}]  (* read a 2x2 matrix in one call *)
In[4]:= Read[ins, Hold[Expression]]  (* read the next expression WITHOUT evaluating it *)
In[5]:= Close[ins]
```

### Notes

`Read` reads **one** object (or one nested type structure) from a stream and
advances its current point, so successive `Read` calls walk the file;
`ReadList` is this primitive looped to end of file. The type may be
`Byte`, `Character`, `Word`, `Record`, `String`, `Number`, `Real`, or
`Expression`, and a *structure* of types (`{{Number, Number}, ...}`, `Hold[...]`,
any head) is filled depth-first — the matrix above is read with a single call.

The `Expression` leaf is read unevaluated and the assembled result is then
evaluated, so `Read[s, Expression]` evaluates while `Read[s, Hold[Expression]]`
keeps the raw form. Past end of file `Read` returns `EndOfFile`; trailing slots of
a partly-read structure are padded with `EndOfFile`. The separator options
(`RecordSeparators`, `WordSeparators`, `TokenWords`, `NullRecords`, `NullWords`)
default from `Options[Read]`, so `SetOptions[Read, ...]` works.
