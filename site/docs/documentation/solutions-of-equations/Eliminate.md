# Eliminate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Eliminate[eqns, vars]`**

eliminates vars between a list/conjunction of simultaneous equations lhs == rhs, returning a balanced Equal\[\] or an And\[\] of Equal\[\]s in the remaining variables (True if the elimination ideal is empty, False if the system is inconsistent). Works on polynomial equations over Q via a lexicographic Gröbner basis with elimination block. A principal-branch inverse-function pre-pass peels single Sin/ Cos/Tan/Sinh/Cosh/Tanh/Exp/Log wrappers and emits Eliminate::ifun; non-polynomial systems otherwise return unevaluated with Eliminate::nlin.

## Examples (20)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (11)

```mathematica
In[1]:= Eliminate[{x == 2 + y, y == z}, y]
Out[1]= 2 + z == x

In[2]:= Eliminate[{f == x^5 + y^5, a == x + y, b == x y}, {x, y}]
Out[2]= a^5 + 5 a b^2 == 5 a^3 b + f

In[3]:= Eliminate[(2 x + 3 y + 4 z == 1) && (9 x + 8 y + 7 z == 2), z]
Out[3]= 22 x + 11 y == 1

In[4]:= Eliminate[{x^2 + y^2 + z^2 == 1, x - y + z == 2, x^3 - y^2 == z + 1}, z] && 12 + 2 x + 5 x^2 + y == 8 x^3 + 4 x^4 + 2 x^5
Out[4]= 27 + 4 x^2 + 8 x^4 + 4 x^5 + 4 x^6 == 18 x + 28 x^3 && 12 + 2 x + 5 x^2 + y == 8 x^3 + 4 x^4 + 2 x^5 && 12 + 2 x + 5 x^2 + y == 8 x^3 + 4 x^4 + 2 x^5

In[5]:= Eliminate[{x^2 + y^2 + z^2 == 1, x - y + z == 2, x^3 - y^2 == z + 1}, {y, z}]
Out[5]= 27 + 4 x^2 + 8 x^4 + 4 x^5 + 4 x^6 == 18 x + 28 x^3
```

Common-root condition

```mathematica
In[6]:= Eliminate[{x - a == 0, x - b == 0}, x]
Out[6]= b == a
```

Inconsistent -> False

```mathematica
In[7]:= Eliminate[1 == 2, x]
Out[7]= False
```

Solvable -> True

```mathematica
In[8]:= Eliminate[x + y == 0, y]
Out[8]= True
```

```mathematica
In[9]:= Eliminate[{u == Sqrt[x^2 + 1], v == 1/Sqrt[x^2 + 1]}, x]
Out[9]= u v == 1
```

Algebraisation pre-pass; equivalent to Mathematica's u^2 Dt[u]^2 + u(-2 Dt[u]^2 - 4 Dt[y]^2) == -Dt[u]^2

```mathematica
In[10]:= Eliminate[{Dt[y] == x^3/Sqrt[x^2 + 1] Dt[x], u == x^2 + 1, Dt[u] == 2 x Dt[x]}, {Dt[x], x}]
Out[10]= Dt[u]^2 + u^2 Dt[u]^2 == 2 u Dt[u]^2 + 4 u Dt[y]^2
```

```mathematica
In[11]:= Eliminate[{u == x^(1/3), v == x^(2/3)}, x]
Out[11]= v == u^2
```

### Scope (3)

```mathematica
In[12]:= HornerForm[11 x^3 - 4 x^2 + 7 x + 2]
Out[12]= 2 + x (7 + x (-4 + 11 x))

In[13]:= HornerForm[a + b x + c x^2, x]
Out[13]= a + x (b + c x)

In[14]:= HornerForm[(11 x^3 - 4 x^2 + 7 x + 2)/(x^2 - 3 x + 1)]
Out[14]= (2 + x (7 + x (-4 + 11 x)))/(1 + x (-3 + x))
```

### Applications (6)

```mathematica
In[15]:= Eliminate[{x + y == 2, x - y == 0}, y]
Out[15]= x == 1

In[16]:= Eliminate[{a == b + c, d == a - c}, c]
Out[16]= d == b

In[17]:= Eliminate[{a == x + y, b == x y}, {x, y}]
Out[17]= True

In[18]:= Eliminate[{p == x + 1/x, q == x^2 + 1/x^2}, x]
Out[18]= 2 + q == p^2

In[19]:= Eliminate[{x == a Cos[t], y == a Sin[t]}, t]
Out[19]= x^2 + y^2 == a^2

In[20]:= Eliminate[{u == Exp[x], v == Exp[2 x]}, x]
Out[20]= v == u^2
```

## Options & behaviour

### Features

- `Protected`.
- Nests multiplications instead of using powers (e.g., $a + x(b + c x)$ instead of $a + bx + cx^2$).
- Identifies variables using `Variables` if not explicitly specified.
- Issues an error and returns unevaluated if the expression is not a polynomial or rational function in the target variables.

## Algorithm

eliminate.c

```text
`Eliminate[eqns, vars]` -- the user-facing variable-elimination front
door.  Takes a list / conjunction of `lhs == rhs` equations and a list
```

of variables to eliminate, then drives the lex-order Buchberger engine in groebner.c with the elimination-block layout that GroebnerBasis's

```text
3-arg form already uses internally.  Surviving basis polynomials are
```

re-presented as balanced `Equal[posPart, negPart]` equations and combined with `And` for multiples.

The pre-pass tries to handle simple transcendental equations of the shape `f[poly] == const` (or `const == f[poly]`) for invertible elementaries f -- one principal-branch rewrite per layer, with an

```text
`Eliminate::ifun` diagnostic to flag that solutions may be missed.

Memory contract: the evaluator owns `res`.  Every early-return path
```

frees all temporaries we own (wrapper Lists, all_vars arrays,

```text
partially-built GBPoly* arrays); we never `expr_free(res)`.  See
```

SPEC.md §4 for the builtin ownership rule.

## Implementation notes

**Algorithm.** `builtin_eliminate` (in `src/poly/eliminate.c`) removes a set of variables from a system of equations by **Gröbner elimination**. It accepts `Eliminate[eqns, vars]` where `eqns` is a `List`/`And` of `lhs == rhs` equations. A pre-pass tries to handle simple invertible transcendental equations of the shape `f[poly] == const` (e.g. `Exp`, `Log`, trig) by a one-layer principal-branch rewrite, emitting an `Eliminate::ifun` diagnostic to warn that branches may be lost.

The algebraic path moves each equation to `lhs − rhs` form, collects all variables, and orders them so the variables to be eliminated form the leading block (the same elimination-block layout used by `GroebnerBasis`'s 3-arg form). It converts the polynomials to `GBPoly` via `gb_from_expr` and runs `gb_buchberger` under the lex/elimination order. By the elimination theorem, the basis polynomials free of the eliminated variables generate the elimination ideal; those survivors are re-presented as balanced `Equal[posPart, negPart]` equations and joined with `And` when there are several.

**Data structures.** Reuses the Gröbner subsystem's `GBPoly` (GMP `mpq_t` coefficients + row-major exponent matrix) and the Buchberger/Gebauer–Möller core in `groebner.c`. The memory contract is the standard builtin one — every early-return path frees temporaries and never frees `res`.

**Complexity / limits.** Inherits Buchberger's worst-case doubly-exponential cost (lex/elimination orders are the expensive case). `gb_from_expr` cannot atomise Power-headed main variables, so genuinely transcendental systems outside the simple-invertible pre-pass fall back to leaving the input unevaluated.

- `Protected`.
- Drives the lex-order Buchberger engine (`GroebnerBasis`) with an
  elimination block: `vars` are placed in the high-priority block so
  that any polynomial mentioning them is filtered out before output.
- A single-layer principal-branch inverse-function pre-pass handles
  equations of the form `f[u] == v` (or `v == f[u]`) for `f` in
  `{Sin, Cos, Tan, Sinh, Cosh, Tanh, Exp, Log}` and their inverses,
  rewriting to `u == InverseF[v]`.  Whenever a rewrite fires Eliminate
  emits the `Eliminate::ifun` diagnostic, signalling that the principal
  branch was used and some solutions may be missed — use `Reduce` for
  complete solution information.
- An algebraisation pre-pass handles inputs containing
  `Power[base, p/q]` with `q > 1` (i.e. `Sqrt[...]`, `(...)^(1/3)`,
  `(...)^(3/2)`, `1/Sqrt[...]`, ...) whose base mentions an elim
  variable.  Each unique base gets a fresh auxiliary symbol `$elN$`
  with the constraint `Power[$elN$, L] == base` (L the LCM of all
  rational-exponent denominators for that base), the power is rewritten
  as `Power[$elN$, p*L/q]`, and `$elN$` is appended to the elim list.
  Nested radicals (e.g. `Sqrt[x + Sqrt[x]]`) work in one shot because
  the constraint substitutes inner aux vars as well.  When the pass
  fires, `Eliminate::alg` is emitted to flag that the returned
  polynomial relation is the cross-multiplied generic consequence —
  sign / branch information may be lost.
- A transcendental algebraisation pre-pass (the general sibling of the
  radical pass) handles elim variables appearing inside circular /
  hyperbolic trig, exponential (`b^x`), or logarithmic functions — even
  in multiple places, where the single-layer inverse pre-pass cannot
  help.  Each equation is first expanded (`TrigExpand` for trig;
  `b^(p+q) -> b^p b^q`, `b^(k m) -> (b^m)^k` for exp; `Log[a b] ->
  Log[a]+Log[b]`, `Log[a^n] -> n Log[a]` for log) so every kernel lands
  on an *atomic* argument; thus `Sin[x]` and `Sin[3x]`, or `Exp[x]` and
  `Exp[2x]`, collapse onto one shared aux.  Exponentials whose exponents
  differ only by a *rational* factor of a shared monomial are commensurate
  and also collapse: for each `(base, monomial)` group the atomic kernel is
  `base^((1/L) m)` with `L` the LCM of the coefficient denominators, so
  `E^(x/3)` and `E^(x/2)` both become powers of `E^(x/6)`
  (`E^(x/3) = (E^(x/6))^2`, `E^(x/2) = (E^(x/6))^3`) — without this they
  register as independent kernels and Gate A rejects the system as `nlin`.  Each `Sin[θ]`/`Cos[θ]`
  (resp. `Sinh`/`Cosh`) becomes a fresh aux pair `$tsN$`,`$tcN$` with
  the Pythagorean constraint `s^2 + c^2 == 1` (circular) /
  `c^2 - s^2 == 1` (hyperbolic); each `b^θ` / `Log[θ]` becomes a single
  algebraically-free aux `$teN$` / `$tlN$`.  The auxes join the elim
  list and `Eliminate::ifun` is emitted (branch / sign information may be
  lost).  Two conservative soundness gates bail to `Eliminate::nlin`
  rather than emit a wrong relation: when an elim variable is shared
  across two distinct kernels (e.g. `Sin[x]` and `Sin[x*y]`, or `Sin[x]`
  and `Exp[x]`), or when it still appears as a bare polynomial atom after
  substitution (e.g. `x` alongside `Sin[x]`).  This is what lets the
  integration-by-substitution pattern
  `Eliminate[{Dt[y] == Cos[x] Sqrt[1-Sin[x]] Dt[x], u == Sin[x],
  Dt[u] == Cos[x] Dt[x]}, {x, Dt[x]}]` return `Dt[y]^2 == (1-u) Dt[u]^2`.
- An inverse-function substitution pre-pass handles the mirror shape,
  where an *inverse* trig / hyperbolic function of an elim variable sits
  buried inside a product (so the single-layer inverse pre-pass cannot
  reach it) but a *defining* equation `M == ArcF[x]` is present (`M`
  elim-free, `x` a single elim symbol, `ArcF` in `{ArcSin, ArcCos,
  ArcTan, ArcSinh, ArcCosh, ArcTanh}`).  The pass propagates
  `ArcF[x] -> M` through the whole system, rewrites the companion base
  `c + s·x^2 = co[M]^(2·co_sign)` to a power of the co-function of the
  *main* angle for **any** rational exponent — the half-integer (radical)
  case `Sqrt[1-x^2] -> Cos[M]` (ArcSin; `-> Sin[M]` ArcCos, `Sqrt[1+x^2]
  -> Cosh[M]` ArcSinh, `Sqrt[x^2-1] -> Sinh[M]` ArcCosh, `-> Sec[M]`/
  `Sech[M]` ArcTan/ArcTanh) **and** the integer (rational-denominator) case
  `1/(1+x^2) -> Cos[M]^2` (ArcTan), `1/(1-x^2) -> Cosh[M]^2` (ArcTanh) — then
  pins `x == F[M]` (F the forward function).  Because the co-function is a
  trig call of the retained variable `M`, any factor of it left after
  elimination is a main-variable monomial that is divided out cleanly (or
  cancels directly in the Groebner step) — so the u-substitution shapes
  `Eliminate[{Dt[y] == x ArcSin[x]/Sqrt[1-x^2] Dt[x], u == ArcSin[x],
  Dt[u] == Dt[x]/Sqrt[1-x^2]}, {x, Dt[x]}]` and
  `Eliminate[{Dt[y] == ArcTan[x]/(1+x^2) Dt[x], u == ArcTan[x],
  Dt[u] == Dt[x]/(1+x^2)}, {x, Dt[x]}]` return `u Sin[u] Dt[u] == Dt[y]`
  and the fully-reduced `u Dt[u] == Dt[y]` respectively (with
  `Eliminate::ifun`).  A `Log` defining
  equation `M == Log[x]` is handled by the same pass (forward function
  `Exp`, no companion radical) — but only as a *fallback*, when `x` also
  occurs as a genuine polynomial atom (e.g. `1/x`), which is exactly the
  case the forward log pass cannot capture (Gate B).  Rather than
  substituting `x -> E^M`, it pins `x == E^M` so `E^M` stays a single
  main-variable atom (the Groebner term reader treats a `Power[base, exp]`
  with symbolic exponent as one opaque indeterminate); when `x` appears
  solely inside logs the forward pass gives a tidier answer and is used
  instead.  So `Eliminate[{Dt[y] == Log[x]/x Dt[x], u == Log[x],
  Dt[u] == -Dt[x]/x^2}, {x, Dt[x]}]` returns `u Dt[u] E^u + Dt[y] == 0`
  (i.e. `-Dt[y] == E^u u Dt[u]`).  This is the pass that unlocks
  `Integrate[x ArcSin[x]/Sqrt[1-x^2], x]` and its inverse-trig /
  inverse-hyperbolic kin — and now `Integrate[Log[x]/x, x]`-style log
  substitutions — through `DerivativeDivides`.
- Equations are normalised through `Numerator[Together[lhs - rhs]]`
  before Buchberger, clearing any `Power[t, -k]` denominators
  introduced by the algebraisation pre-passes.  Surviving basis
  polynomials are stripped of any common monomial factor so the
  reported equation matches the primitive form (no spurious factor
  of a main variable from a `u == Sqrt[...]` equation).
- Equations that remain non-polynomial after both pre-passes (e.g.
  `Sin[x*y] == ...` with `x` in the elim set) return the expression
  unevaluated with `Eliminate::nlin`.
- Main-block variables (the implicit free parameters) are discovered
  automatically by walking the equations; mathematical constants such
  as `Pi`, `E`, and `EulerGamma` flow through as parameter symbols.
  Function-shaped subexpressions whose head is not one of
  `{Plus, Times, Power, List, And, Equal, Or, Not}` are treated as
  single polynomial atoms (so `Dt[y]`, `f[a]`, `Sin[k]` etc. are
  whole-expression variables rather than products of their heads and
  arguments).

**Attributes:** `Protected`.

## References

**See also:** [GroebnerBasis](../../algebra/GroebnerBasis/), [Reduce](../../solutions-of-equations/Reduce/), [TrigExpand](../../elementary-functions/TrigExpand/), [Sinh](../../elementary-functions/Sinh/), [Cosh](../../elementary-functions/Cosh/), [Log](../../elementary-functions/Log/), [Exp](../../elementary-functions/Exp/), [Pi](../../mathematical-constants/Pi/)

- T. Becker, V. Weispfenning, *Gröbner Bases* (Springer, 1993).
- D. Cox, J. Little, D. O'Shea, *Ideals, Varieties, and Algorithms* (Springer).
- Source: [`src/poly/eliminate.c`](https://github.com/stblake/mathilda/blob/main/src/poly/eliminate.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_eliminate.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eliminate.c)

## Notes & additional examples

### Notes

`Eliminate[eqns, vars]` removes `vars` from a system of polynomial equations
over the rationals (via a lexicographic Gröbner basis with an elimination
block), returning the relations among the remaining variables. It yields
`True` if the elimination ideal is empty and `False` if the system is
inconsistent; non-polynomial systems return unevaluated with `Eliminate::nlin`.
