# Real

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Real`**

is the head of approximate real numbers. As a type specification in Read and ReadList it reads a number, always returned as an approximate number (C/Fortran E notation is accepted).

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

Real is the head of approximate numbers

```mathematica
In[1]:= Head[3.14]
Out[1]= Real
```

As a read type: E-notation accepted

```mathematica
In[2]:= s = OpenWrite["/tmp/mathilda_real.txt"]; WriteString[s, "6.022e23"]; Close[s]; Read["/tmp/mathilda_real.txt", Real]
Out[2]= 6.022e+23
```

## Implementation notes

**Definition.** `Real` has two inert roles, both without a builtin. First, it is the
**head of approximate real numbers**: `Head[3.14]` is `Real`, so `Real` names the type
of machine- and arbitrary-precision floating-point values (matched by `_Real` in
patterns). Second, it is a **type specification** for `Read` and `ReadList`, telling the
reader to read a number and always return it as an approximate (inexact) number, even
when the token has no decimal point.

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_Real`). As a value head it is
attached to every `EXPR_REAL` / `EXPR_MPFR` node (the printer shows these as `3.14`,
not `Real[...]`). As a read type, `read.c`'s `type_from_symbol` maps it to the internal
`RT_REAL` reader, which accepts C/Fortran `E`-notation. The related `Number` type reads
an integer when the token is integral and a real otherwise; `Real` always coerces to a
real.

**Usage & limits.** As a read type it is meaningful only inside a `Read`/`ReadList`
call; on its own `Real` evaluates to itself. Past end of file the reader returns
`EndOfFile`.

**Attributes:** none registered.

## References

- Source: [`src/io/read.c`](https://github.com/stblake/mathilda/blob/main/src/io/read.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)

## Notes & additional examples

### Notes

`Real` has two roles. It is the **head of approximate real numbers**, so `Head[3.14]` is
`Real` and the pattern `_Real` matches any machine- or arbitrary-precision float. It is
also a **type specification** for `Read` and `ReadList`, where it reads a number and
always returns it as an approximate number (accepting C/Fortran `E`-notation), even when
the token is integral.

The related `Number` read type returns an integer when the token has no decimal point or
exponent and an approximate number otherwise; `Real` always coerces to a real. As a read
type it is meaningful only inside a `Read`/`ReadList` call; on its own `Real` evaluates
to itself.
