### Worked examples

```mathematica
In[1]:= Head[Word]  (* a bare read-type token, not a function *)
```

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_word.txt"]; WriteString[str, "alpha beta gamma\n"]; Close[str];
In[2]:= ReadList["/tmp/mathilda_word.txt", Word]  (* each whitespace-delimited run as a string *)
```

```mathematica
In[1]:= ins = OpenRead["/tmp/mathilda_word.txt"]; w = Read[ins, Word]; Close[ins]; w  (* just the first word *)
```

### Notes

`Word` is a **read type** for `Read` and `ReadList`: it reads the next run of characters
delimited by word separators and returns it as a string. It is an inert token symbol — no
builtin, no DownValues — recognised by the reader in `src/io/read.c`, which maps it to its
internal word case.

`ReadList["file", Word]` returns every word as a flat list of strings; `Read[stream, Word]`
returns just the next one and advances the stream. The separators default to space and tab
and are set by `WordSeparators`; `TokenWords` forces given strings to read as standalone
words, and `NullWords -> True` keeps empty fields between adjacent separators. A `Word` read
past end of file returns `EndOfFile`.
