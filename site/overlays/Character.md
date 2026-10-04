### Worked examples

```mathematica
In[1]:= s = OpenWrite["/tmp/mathilda_char.txt"]; WriteString[s, "Hi"]; Close[s]; Read["/tmp/mathilda_char.txt", Character]  (* reads one character as a string *)
```

```mathematica
In[1]:= MemberQ[{Byte, Character, Expression, Number, Real, Record, String, Word}, Character]  (* one of the eight read types *)
```

### Notes

`Character` is a type specification for `Read` and `ReadList`: it consumes exactly one
character from the stream and returns it as a one-character string. It has no value of
its own, so on its own it just evaluates to itself.

It may stand alone or nest inside a type structure (`{Character, Character}`,
`Hold[...]`, any head), which the reader fills depth-first. Contrast `Byte`, which
returns the raw integer code 0-255 instead of a string. Past end of file the reader
returns `EndOfFile`.
