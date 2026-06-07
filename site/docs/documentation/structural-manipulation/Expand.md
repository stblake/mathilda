# Expand

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Expand[expr] expands out products and powers in expr.
Expand[expr, patt] leaves unexpanded any parts of expr that are free of the pattern patt.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Expand[(x+3)(x+2)]
Out[1]= 6 + 5 x + x^2

In[2]:= Expand[(x+y)^2 (x-y)^2]
Out[2]= x^4 - 2 x^2 y^2 + y^4

In[3]:= Expand[(x+1)^2 + (y+1)^2, x]
Out[3]= 1 + 2 x + x^2 + (1 + y)^2
```

## Implementation notes

**Algorithm.** `builtin_expand` (in `src/expand.c`) calls `expr_expand_patt(e, patt)`, a structural recursion that distributes products over sums. The dispatch by head:

- **Plus** — expand each summand and rebuild via `eval_and_free`.
- **Times** — expand each factor, then multiply them pairwise with `multiply_all`, a divide-and-conquer fold whose leaf operation `multiply_two` handles the four `Plus×Plus / Plus×atom / atom×Plus / atom×atom` cases, producing every cross term (`a_i · b_j`) and re-summing them.
- **Power[base, k]** with `k` a positive integer `< 100` — expand the base, then `power_expand` raises it by **repeated squaring** (binary exponentiation), each multiply going through the distributing `multiply_two`. (Note: there is no explicit binomial-coefficient formula; the binomial expansion of `(a+b)^k` falls out of the iterated distribution.) Negative or non-integer exponents are left untouched.
- **List / equations / inequalities / And / Or / Not** — threads into each argument (passing operator-symbol slots of `Inequality` through verbatim).

A two-argument `Expand[expr, patt]` only expands subexpressions that *contain* `patt` (gated by `expr_contains_patt`, which uses the pattern matcher); everything else is copied unchanged. Expand does **not** descend into the arguments of arbitrary function heads.

**Data structures.** Plain `Expr*` trees; transient `Expr**` argument buffers per node. All arithmetic rebuilding goes back through the evaluator (`eval_and_free`) so `Plus`/`Times` canonicalisation (Flat/Orderless, like-term collection) applies after each step.

**Complexity / limits.** Worst case is the full cross-product blow-up of distribution: expanding `(x_1+...+x_m)^k` produces `O(m^k)`-sized output, and `Power` is capped at exponent `< 100` to bound it.

- `Protected`.
- Works only on positive integer powers and distributes products over sums.
- Threads over equations, inequalities, and lists.
- Implements an efficient binary-splitting algorithm for distributing products and repeated squaring for powers.
- `Expand[expr, patt]` leaves unexpanded any parts of `expr` that are free of the pattern `patt`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), Ch. 3 (normal forms and the distributive expansion of polynomials).
- Source: [`src/expand.c`](https://github.com/stblake/mathilda/blob/main/src/expand.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Expand[(x + 1)^3]
Out[1]= 1 + 3 x + 3 x^2 + x^3
```

```mathematica
In[1]:= Expand[(x + 1)^4]
Out[1]= 1 + 4 x + 6 x^2 + 4 x^3 + x^4
```

```mathematica
In[1]:= Expand[(a + b)(c + d)]
Out[1]= a c + b c + a d + b d
```

```mathematica
In[1]:= Expand[(x + 2)^2 (x - 1)]
Out[1]= -4 + 3 x^2 + x^3
```

### Notes

`Expand` applies the distributive law to products and integer powers,
producing a flat sum of monomials in canonical order (ascending total
degree in the leading variable). Like terms are combined automatically, so
`(x + 2)^2 (x - 1)` collapses the `x^1` coefficient to zero and it drops out
of the result. `Expand` only multiplies out — it does not factor or cancel —
and a second argument `Expand[expr, patt]` leaves alone any parts free of the
pattern `patt`.
