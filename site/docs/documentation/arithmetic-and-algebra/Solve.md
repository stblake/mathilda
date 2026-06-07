# Solve

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Solve[expr, vars]
    Attempts to solve the equation or system expr for the
    variables vars.
Solve[expr, vars, dom]
    Solves over the domain dom.  Default Complexes; Reals filters
    down to real roots via per-degree discriminant and sign tests;
    Integers further restricts the output to provably concrete
    integer solutions (Integer / BigInt only -- Rationals, Sqrt[],
    and held Root[] objects are dropped).

Options:
    Cubics              -> False     (radical form for cubics)
    Quartics            -> False     (radical form for quartics)
    InverseFunctions    -> Automatic (use inverse-function peel)
    GeneratedParameters -> C         (head for parameters C[k])
    VerifySolutions     -> Automatic (reserved)

Solves single polynomial equalities, radical equations, linear
systems, and -- via the inverse-function specialist -- single-
variable equations whose outermost dependence is an elementary
invertible head (Log, Exp, Sin/Cos/Tan/Cot/Sec/Csc, their
hyperbolic counterparts, the inverse trig/hyperbolic forms,
and Power[g, n] for integer n >= 2).  Multi-branch heads
introduce an integer parameter C[k] wrapped in
ConditionalExpression[..., Element[C[k], Integers]].  Emits
Solve::ifun the first time inverse functions are used.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Solve[2 x + 3 == 0, x]
Out[1]= {{x -> -3/2}}

In[2]:= Solve[x^2 - 5 x + 6 == 0, x]
Out[2]= {{x -> 2}, {x -> 3}}

In[3]:= Solve[x^2 + 1 == 0, x]
Out[3]= {{x -> -I}, {x -> I}}

In[4]:= Solve[x^2 + 1 == 0, x, Reals]
Out[4]= {}

In[5]:= Solve[(x-1)^2 == 0, x]
Out[5]= {{x -> 1}, {x -> 1}}

In[6]:= Solve[x^4 - 5 x^2 + 4 == 0, x]
Out[6]= {{x -> -2}, {x -> -1}, {x -> 1}, {x -> 2}}

In[7]:= Solve[x^3 + x + 1 == 0, x]
Out[7]= {{x -> Root[1 + #1 + #1^3 &, 1]}, {x -> Root[1 + #1 + #1^3 &, 2]}, {x -> Root[1 + #1 + #1^3 &, 3]}}

In[8]:= Solve[Sin[x] == 0, x]
Out[8]= {{x -> ConditionalExpression[Pi + 2 C[1] Pi, Element[C[1], Integers]]}, {x -> ConditionalExpression[2 C[1] Pi, Element[C[1], Integers]]}}
```

## Implementation notes

- `Protected`.  Matches Mathematica's attribute set -- arguments are

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- von zur Gathen & Gerhard, "Modern Computer Algebra" (3rd ed.), Ch. 14 (polynomial roots and resolution).
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), Ch. 9 (solving systems).
- Source: [`src/solve.c`](https://github.com/stblake/mathilda/blob/main/src/solve.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Solve[x^2 - 5 x + 6 == 0, x]
Out[1]= {{x -> 2}, {x -> 3}}
```

```mathematica
In[1]:= Solve[x^2 + 1 == 0, x]
Out[1]= {{x -> -I}, {x -> I}}
```

```mathematica
In[1]:= Solve[x^2 - 2 == 0, x]
Out[1]= {{x -> -Sqrt[2]}, {x -> Sqrt[2]}}
```

```mathematica
In[1]:= Solve[{x + y == 3, x - y == 1}, {x, y}]
Out[1]= {{x -> 2, y -> 1}}
```

### Notes

`Solve` returns a list of solution rule-lists, one per solution; each inner
list assigns every requested variable. Complex roots are produced by default,
so `x^2 + 1 == 0` yields the conjugate pair `±I`, and irrational roots come
back in exact radical form (`±Sqrt[2]`). Linear systems are solved directly
and return a single rule-list. Cubic roots are reported using `(-1)^(1/3)`
style radicals; pass `Cubics -> False` / `Quartics -> False` to suppress
explicit radical forms when desired.
