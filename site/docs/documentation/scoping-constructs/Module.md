# Module

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Module[{x, y, ...}, expr] specifies that x, y, ... are local variables.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= x = 1; Module[{x = 2}, x + 1]
Out[1]= 3

In[2]:= x
Out[2]= 1
```

## Implementation notes

- `HoldAll`, `Protected`.
- Variables are renamed to `name$nnn` using `$ModuleNumber`.
- Created symbols have the `Temporary` attribute.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Harold Abelson and Gerald Jay Sussman, *Structure and Interpretation of Computer Programs*, 2nd ed., §3.1 (local state and lexical scoping).
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/scoping-constructs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/scoping-constructs.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Module[{x = 5}, x^2 + 1]
Out[1]= 26
```

```mathematica
In[1]:= Module[{a = 2, b = 3}, a*b + a + b]
Out[1]= 11
```

```mathematica
In[1]:= f[n_] := Module[{s = 0}, s = n^2 + n; s]; f[4]
Out[1]= 20
```

### Notes

`Module[{vars}, body]` introduces lexically scoped local variables, optionally
with initial values (`{x = 5}`). The locals are renamed to unique symbols so they
never collide with global definitions of the same name. Inside the body, locals
can be reassigned (`s = n^2 + n`) and a sequence of statements separated by `;`
is evaluated left to right, with the last expression returned. This makes
`Module` the standard tool for writing multi-step function definitions.
