# CylindricalDecomposition

!!! warning "Status: Partial"
    implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## Description

**`CylindricalDecomposition[expr, vars]`**

Gives a cylindrical algebraic decomposition of the real solution set of expr -- a logical combination of polynomial equations and inequalities -- as a quantifier-free And/Or formula in which each variable is bounded cylindrically in terms of the earlier ones, e.g. CylindricalDecomposition\[x^2 + y^2 \<= 1, {x, y}\] gives -1 \<= x \<= 1 && -Sqrt\[1 - x^2\] \<= y \<= Sqrt\[1 - x^2\].  The domain is always the Reals.  Returns True / False when the statement decides, and stays unevaluated when the decomposition cannot be computed exactly (an undecidable sign, or a positive-dimensional system with irrational fibres).  Reduce's options (Modulus, Cubics, Quartics, WorkingPrecision, ...) may be given and are forwarded.

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Worked examples (1)

```mathematica
In[1]:= CylindricalDecomposition[x^2 == -1, {x}]
Out[1]= False
```

### Applications (7)

```mathematica
In[2]:= CylindricalDecomposition[x^2 > 1, {x}]
Out[2]= x < -1 || x > 1

In[3]:= CylindricalDecomposition[x^3 - x > 0, {x}]
Out[3]= -1 < x < 0 || x > 1

In[4]:= CylindricalDecomposition[x y > 1, {x, y}]
Out[4]= x < 0 && y < 1/x || x > 0 && y > 1/x

In[5]:= CylindricalDecomposition[x^2 + y^2 <= 1, {x, y}]
Out[5]= -1 <= x <= 1 && -1/2 Sqrt[4 - 4 x^2] <= y <= 1/2 Sqrt[4 - 4 x^2]

In[6]:= CylindricalDecomposition[x^2 + 1 > 0, {x}]
Out[6]= True

In[7]:= CylindricalDecomposition[x^2 < 0, {x}]
Out[7]= False

In[8]:= CylindricalDecomposition[x^2 == 1, {x}]
Out[8]= x == -1 || x == 1
```

## Implementation notes

**Algorithm.** `builtin_cylindrical_decomposition` is a thin, Reals-only companion to
`Reduce`: it does no decomposition of its own but forwards the problem to `Reduce`'s
real-domain engine, where the cylindrical algebraic decomposition actually lives. This
keeps a single implementation of the CAD/real-quantifier-elimination machinery and
guarantees the two heads agree on every real problem.

*Argument handling.* It first peels any trailing option `Rule`s (a symbol LHS, e.g.
`Modulus -> p`) off the end, leaving the positional arguments. It accepts `[expr, vars]`
or `[expr, vars, Reals]` -- since the domain is always the reals, an explicit `Reals` is
redundant but allowed, and any other third positional declines soundly (returns `NULL`,
leaving the input unevaluated). The variable specification is validated by
`cad_valid_vars`.

*Delegation.* It builds `Reduce[expr, vars, Reals, <peeled options...>]`, forwarding the
option rules verbatim so that all of `Reduce`'s options (`Modulus`, `Cubics`, `Quartics`,
`WorkingPrecision`, ...) are reused with no per-option logic here, and evaluates it under
`mth_msg_suppress_push` / `pop` so that any diagnostics emitted while `Reduce` probes
internally are not misattributed to `CylindricalDecomposition`.

*Sound decline.* If `Reduce` declines -- returning an unevaluated `Reduce[...]` -- the
companion returns `NULL` rather than echoing an inner `Reduce[...]` under its own head,
so `CylindricalDecomposition[...]` is simply left unevaluated. Otherwise `Reduce`'s
result (a quantifier-free `And`/`Or` formula, or a bare `True` / `False`) is returned
directly.

**Complexity / limits.** Entirely those of the delegate `Reduce` over the reals: CAD is
doubly exponential in the number of variables in the worst case, so the head is practical
for few variables and low-degree polynomials, and returns unevaluated where the real
decomposition cannot be produced exactly.

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [Reduce](../../solutions-of-equations/Reduce/), [Modulus](../../other-advanced/Modulus/), [Cubics](../../solutions-of-equations/Cubics/), [Quartics](../../solutions-of-equations/Quartics/)

- Collins, "Quantifier Elimination for Real Closed Fields by Cylindrical Algebraic Decomposition", LNCS 33 (1975).
- Collins & Hong, "Partial Cylindrical Algebraic Decomposition for Quantifier Elimination", J. Symbolic Computation 12 (1991).
- McCallum, "An Improved Projection Operation for Cylindrical Algebraic Decomposition", in Caviness & Johnson (eds.), Quantifier Elimination and CAD (1998).
- Basu, Pollack & Roy, "Algorithms in Real Algebraic Geometry" (2nd ed., 2006), Ch. 5 & 11.
- Source: [`src/solve/reduce_companions.c`](https://github.com/stblake/mathilda/blob/main/src/solve/reduce_companions.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)

## Notes & additional examples

### Notes

`CylindricalDecomposition[expr, vars]` returns a *cylindrical algebraic decomposition*
of the real solution set of `expr` -- a quantifier-free `And`/`Or` formula in which the
variables are bounded one after another, each in terms of the earlier ones, so that the
region is described as a stack of cylinders. It is the geometric companion of
`Reduce` over the reals: the two return the same formula for the same real problem, and
`CylindricalDecomposition` exists to name the decomposition directly and to keep the
domain fixed at the `Reals` (an explicit `Reals` third argument is accepted but
redundant).

Because the domain is always the reals, the endpoints are real algebraic numbers and may
be irrational surds, as in the disk example above. A statement that decides returns
`True` or `False`; one whose decomposition cannot be computed exactly -- an undecidable
sign, or a positive-dimensional system with irrational fibres beyond the current engine
-- is returned unevaluated rather than guessed, so a formula that *is* returned describes
the whole real solution set exactly. All of `Reduce`'s options (`Modulus`, `Cubics`,
`Quartics`, `WorkingPrecision`, ...) may be given and are forwarded unchanged.
