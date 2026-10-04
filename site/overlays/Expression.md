### Worked examples

```mathematica
In[1]:= s = OpenWrite["/tmp/mathilda_expr.txt"]; WriteString[s, "1 + 2*3"]; Close[s]; Read["/tmp/mathilda_expr.txt", Expression]  (* reads and evaluates one expression *)
```

```mathematica
In[1]:= MemberQ[{Byte, Character, Expression, Number, Real, Record, String, Word}, Expression]  (* one of the eight read types *)
```

### Notes

`Expression` is a type specification for `Read` and `ReadList`: it parses one complete
Mathilda expression from the stream. It is the **default** read type — with no type
given, `ReadList["file"]` reads every remaining expression.

The leaf is read unevaluated and the assembled result is then evaluated, so
`Read[s, Expression]` evaluates what it read (the example above returns `7`) while
`Read[s, Hold[Expression]]` keeps the raw, unevaluated form. On its own `Expression`
evaluates to itself; past end of file the reader returns `EndOfFile`.
