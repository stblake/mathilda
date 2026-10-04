# Expression

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Expression`**

is a type specification in Read and ReadList that reads one complete Mathilda expression.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

Reads and evaluates one expression

```mathematica
In[1]:= s = OpenWrite["/tmp/mathilda_expr.txt"]; WriteString[s, "1 + 2*3"]; Close[s]; Read["/tmp/mathilda_expr.txt", Expression]
Out[1]= 7
```

One of the eight read types

```mathematica
In[2]:= MemberQ[{Byte, Character, Expression, Number, Real, Record, String, Word}, Expression]
Out[2]= True
```

## Implementation notes

**Definition.** `Expression` is a **type specification** used by `Read` and `ReadList`.
It has no builtin and no value — it is an inert token that tells the reader to parse one
complete Mathilda expression from the stream. It is the **default** type when none is
given (both `read.c` and `readlist.c` fall back to `SYM_Expression`).

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_Expression`). `read.c`'s
`type_from_symbol` maps it to the internal `RT_EXPRESSION` reader. The leaf is read
unevaluated and the assembled result is then evaluated, so `Read[s, Expression]`
evaluates the expression it read while `Read[s, Hold[Expression]]` keeps it in raw form.

**Usage & limits.** Meaningful only as a `Read`/`ReadList` type argument; evaluated on
its own it returns itself. With no explicit type, `ReadList["file"]` reads every
remaining expression. Past end of file the reader returns `EndOfFile`.

**Attributes:** none registered.

## References

- Source: [`src/io/read.c`](https://github.com/stblake/mathilda/blob/main/src/io/read.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`Expression` is a type specification for `Read` and `ReadList`: it parses one complete
Mathilda expression from the stream. It is the **default** read type — with no type
given, `ReadList["file"]` reads every remaining expression.

The leaf is read unevaluated and the assembled result is then evaluated, so
`Read[s, Expression]` evaluates what it read (the example above returns `7`) while
`Read[s, Hold[Expression]]` keeps the raw, unevaluated form. On its own `Expression`
evaluates to itself; past end of file the reader returns `EndOfFile`.
