### Worked examples

```mathematica
In[1]:= s = OpenWrite["/tmp/mathilda_tok.txt"]; WriteString[s, "a+b+c"]; Close[s]; ReadList["/tmp/mathilda_tok.txt", Word, TokenWords -> {"+"}]  (* "+" read as its own word *)
```

### Notes

`TokenWords` is an option for `Read` and `ReadList` giving strings that are read as
separate words even when they are not surrounded by word separators. The default is `{}`.

With `TokenWords -> {"+"}`, the `Word` reader emits `+` as its own field, so `"a+b+c"`
reads as `{"a", "+", "b", "+", "c"}` — the standard way to tokenise an expression whose
operators are glued to their operands. It is an inert, `Protected` keyword with a
registered default, so `SetOptions[Read, TokenWords -> ...]` works. Related tokenisation
options are `WordSeparators`, `RecordSeparators`, `NullWords` and `NullRecords`.
