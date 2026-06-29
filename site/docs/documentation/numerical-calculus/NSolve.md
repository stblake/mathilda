# NSolve

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NSolve[expr, vars]
    gives numerical approximations to the solutions of the equation or system expr for the variables vars, as a list of replacement-rule lists. NSolve[expr, vars, Reals] restricts to real solutions; the default domain is the complexes. vars may be a single variable or a list; NSolve[{e1, e2, ...}, vars] is the conjunction e1 && e2 && .... A working precision may be given as a trailing positional argument or via WorkingPrecision. Results: {} no solutions, {{x->s,...},...} the solutions (univariate roots are repeated by multiplicity), {{}} the universal solution. A univariate polynomial equation is solved with NRoots; square zero-dimensional polynomial systems use a Groebner-basis multiplication-matrix eigenvalue method (Method -> "Symbolic" uses lexicographic elimination); other equations fall back to Solve or FindRoot seeding. Integer, real, and complex coefficients are handled at machine and arbitrary precision.

Options: MaxRoots, Method (Automatic | "EndomorphismMatrix" | "Homotopy" | "Symbolic"), WorkingPrecision, VerifySolutions, RandomSeeding.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= NSolve[x^5 - 2 x + 3 == 0, x, Reals]
Out[1]= {{x -> -1.42361}}

In[2]:= NSolve[{x^2 + y^2 == 1, x^3 - y^3 == 2}, {x, y}]
Out[2]= {{x -> -1.09791 + 0.839887*I, y -> 1.09791 + 0.839887*I}, {x -> -1.09791 - 0.839887*I, y -> 1.09791 - 0.839887*I}, {x -> 1.22333 - 0.0729987*I, y -> 0.125423 + 0.712005*I}, {x -> 1.22333 + 0.0729987*I, y -> 0.125423 - 0.712005*I}, {x -> -0.125423 + 0.712005*I, y -> -1.22333 - 0.0729987*I}, {x -> -0.125423 - 0.712005*I, y -> -1.22333 + 0.0729987*I}}

In[3]:= NSolve[{x^2 + y^3 == 1, 2 x + 3 y == 4}, {x, y}, Reals]
Out[3]= {{x -> 7.93641, y -> -3.95761}}

In[4]:= NSolve[x + 2 y + 3 z == 4 && 3 x + 4 y + 5 z == 6 && 6 x + 7 y + 8 z == 0, {x, y, z}]
Out[4]= {}

In[5]:= NSolve[E^x - x == 7, x, Reals]
Out[5]= {{x -> -6.99909}, {x -> 2.22154}}

In[6]:= NSolve[{x^2 + y^2 == 1, x^3 - y^3 == 2}, {x, y}, WorkingPrecision -> 25]
Out[6]= {{x -> -1.0979116727228235764163996 + 0.83988692161565920362280281*I, y -> 1.0979116727228235764163996 + 0.83988692161565920362280281*I}, {x -> -1.0979116727228235764163996 - 0.83988692161565920362280281*I, y -> 1.0979116727228235764163996 - 0.83988692161565920362280281*I}, {x -> 1.2233348984131033766895813 - 0.072998738390442569855466144*I, y -> 0.12542322569027980027318178 + 0.71200452485314764855498901*I}, {x -> 1.2233348984131033766895813 + 0.072998738390442569855466144*I, y -> 0.12542322569027980027318178 - 0.71200452485314764855498901*I}, {x -> -0.12542322569027980027318178 + 0.71200452485314764855498901*I, y -> -1.2233348984131033766895813 - 0.072998738390442569855466144*I}, {x -> -0.12542322569027980027318178 - 0.71200452485314764855498901*I, y -> -1.2233348984131033766895813 + 0.072998738390442569855466144*I}}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
