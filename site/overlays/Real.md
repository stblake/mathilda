### Worked examples

```mathematica
In[1]:= Head[3.14]  (* Real is the head of approximate numbers *)
```

```mathematica
In[1]:= s = OpenWrite["/tmp/mathilda_real.txt"]; WriteString[s, "6.022e23"]; Close[s]; Read["/tmp/mathilda_real.txt", Real]  (* as a read type: E-notation accepted *)
```

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
