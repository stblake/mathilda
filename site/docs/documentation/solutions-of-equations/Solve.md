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
systems, zero-dimensional nonlinear polynomial systems (via a
lexicographic Groebner basis and triangular back-substitution;
positive-dimensional systems emit Solve::nsdim and stay
unevaluated), and -- via the inverse-function specialist --
single-variable equations whose outermost dependence is an elementary
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

**Algorithm.** `builtin_solve` is a classifier/router, not a solver: it parses options, validates the variable spec, normalises the input, then dispatches to one of five specialists in `src/solvepoly.c`, `src/solverad.c`, `src/solvelinsys.c`, `src/solvetrig.c`, `src/solveinv.c`.

The router first peels trailing `Rule`/`RuleDelayed` options (`Cubics`, `Quartics`, `InverseFunctions`, `GeneratedParameters`, `VerifySolutions`, `Assumptions`, `Method`, `Modulus`) off the end of the positional args via `is_known_option_name`; an unrecognised trailing option name emits `Solve::optx`. Positional args are `expr [, vars [, dom]]`. `is_valid_solve_vars` rejects numeric-literal variables with `Solve::ivar`. Compound variables (`Dt[y]`, `f[a,b]`, `x^2`) are rewritten to fresh internal symbols `Solve$var$N` throughout `expr` by `collect_and_subst_compound_vars`, then restored in the output by `unsubst_compound_vars`, so the specialists only ever see bare symbols. Inexact-coefficient inputs are detected (`common_scan_inexact`), force-rationalised to the minimum bit precision found (`common_rationalize_input`), solved exactly, and numericalised back at the tail (`common_numericalize_result`) — exact-in/exact-out, inexact-in/inexact-out, mirroring Cancel/Together/Integrate. `True`/`False` short-circuit to `{{}}` (tautology) / `{}` (contradiction). `Abs[u]==0` is rewritten to `u==0` (`try_abs_zero_rewrite`).

The dispatch cascade: an `And`/`List` of equations, or a single `Equal` over a ≥2-symbol variable list, routes to `solvelinsys_solve_linear_system` (linear-system specialist, canonicalises each equation to `lhs - rhs`, returns NULL when non-affine). Otherwise the single-variable path tries `solvepoly_solve_polynomial_equality` first (the polynomial specialist, `src/poly/solvepoly.c`, also exposed as `Solve\`SolvePolynomialEquality`); on NULL it falls back in order to `solveinv_solve_inverse_equality` (peels one elementary invertible head — Log/Exp/trig/hyperbolic/inverse-trig and integer Power — introducing `C[k]` integer parameters wrapped in `ConditionalExpression`), then `solvetrig_solve_trig_equality` (multi-trig canonicalisation), then `solverad_solve_radicals_equality` (radical equations). A specialist returning NULL leaves the call unevaluated.

The `dom` third argument selects the solution domain: default `Complexes`; `Reals` filters via per-degree discriminant/sign tests inside the polynomial specialist; `Integers` further drops non-concrete-integer solutions. `Cubics`/`Quartics -> False` (the defaults) return cubic/quartic roots as held `Root[]` objects rather than radical formulas.

**Data structures.** Everything is `Expr*`. Options accumulate into a stack `SolveOpts` bundling `SolvePolyOpts` and `SolveInvOpts`. Compound-var substitutions are tracked in a fixed `SolveVarSub subs[32]` array (`SOLVE_MAX_VAR_SUBS`). Output is the standard `List` of `List` of `Rule` rewrite-rule form `{{x -> ...}, ...}`.

**Complexity / limits.** Dominated by the chosen specialist (polynomial root-finding, linear-system elimination, radical isolation). The router itself is linear in expression size plus the substitution passes. The compound-variable cap is 32 distinct variables per call.

- `Protected`.  Uses the standard attribute set -- arguments are

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- von zur Gathen & Gerhard, "Modern Computer Algebra" (3rd ed.), Ch. 14 (polynomial roots and resolution).
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), Ch. 9 (solving systems).
- Source: [`src/solve.c`](https://github.com/stblake/mathilda/blob/main/src/solve.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)

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

```mathematica
In[1]:= Solve[x y == 1 && x + y == 3, {x, y}]
Out[1]= {{x -> 1/2 (3 - Sqrt[5]), y -> 1/2 (3 + Sqrt[5])}, {x -> 1/2 (3 + Sqrt[5]), y -> 1/2 (3 - Sqrt[5])}}

In[2]:= Solve[x y == 6 && x + y == 5, {x, y}, Integers]
Out[2]= {{x -> 3, y -> 2}, {x -> 2, y -> 3}}
```

```mathematica
In[1]:= Solve[a x^2 + b x + c == 0, x]
Out[1]= {{x -> (1/2 (-b + Sqrt[b^2 - 4 a c]))/a}, {x -> (1/2 (-b - Sqrt[b^2 - 4 a c]))/a}}
```

```mathematica
In[1]:= Solve[x^4 - 1 == 0, x]
Out[1]= {{x -> -1}, {x -> 1}, {x -> -I}, {x -> I}}

In[2]:= Solve[Sin[x] == 0, x]
Out[2]= {{x -> ConditionalExpression[Pi + 2 C[1] Pi, Element[C[1], Integers]]}, {x -> ConditionalExpression[2 C[1] Pi, Element[C[1], Integers]]}}
```

### Notes

`Solve` returns a list of solution rule-lists, one per solution; each inner
list assigns every requested variable. Complex roots are produced by default,
so `x^2 + 1 == 0` yields the conjugate pair `±I`, and irrational roots come
back in exact radical form (`±Sqrt[2]`). Linear systems are solved directly
and return a single rule-list. Cubic roots are reported using `(-1)^(1/3)`
style radicals; pass `Cubics -> False` / `Quartics -> False` to suppress
explicit radical forms when desired.

Nonlinear polynomial systems with finitely many solutions (a zero-dimensional
ideal) are solved via a lexicographic Gröbner basis and triangular
back-substitution, honouring the `Reals` / `Integers` domain. Systems with
infinitely many solutions (positive-dimensional ideals, e.g.
`Solve[x^2 - y^2 == 0, {x, y}]`) emit `Solve::nsdim` and are left unevaluated.
