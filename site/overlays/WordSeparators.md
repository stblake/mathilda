### Worked examples

```mathematica
In[1]:= WordSeparators -> {" ", ","}  (* an inert option keyword: it just names the left of a rule *)
```

```mathematica
In[1]:= Options[ReadList]  (* its default is {" ", "\t"}, alongside the other read options *)
```

```mathematica
In[1]:= Attributes[WordSeparators]  (* a plain inert symbol, with no value of its own *)
```

### Notes

`WordSeparators` is an option for `Read` and `ReadList` giving the strings that
separate words when the type specification is `Word`; the default is `{" ", "\t"}`.
It is an inert keyword — the reader consults its value when cutting `Word` tokens
from a stream. A read that uses it, such as `ReadList["data.txt", Word,
WordSeparators -> {",", " "}]`, needs a file (or open stream) on disk, so the
runnable examples above show it only in isolation — as the left of an option rule
and as the default surfaced by `Options[ReadList]`. It does not affect `Number`,
`Record` (which uses `RecordSeparators`), or whole-line `String` reads.
