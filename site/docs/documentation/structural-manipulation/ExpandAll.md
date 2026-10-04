# ExpandAll

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ExpandAll[expr]`**

expands out all products and integer powers in any part of expr.

**`ExpandAll[expr, patt]`**

avoids expanding parts of expr that do not contain terms matching patt.

<details>
<summary>Notes</summary>

ExpandAll effectively maps Expand and ExpandDenominator onto every part of expr, including function heads, arguments, exponents, and denominators. ExpandAll automatically threads over lists, as well as equations, inequalities, and logic functions.

</details>

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= ExpandAll[1/(1+x)^3 + Sin[(1+x)^3]]
Out[1]= 1/(1 + 3 x + 3 x^2 + x^3) + Sin[1 + 3 x + 3 x^2 + x^3]

In[2]:= ExpandAll[(x+z)^2/(x+y)^2]
Out[2]= x^2/(x^2 + 2 x y + y^2) + 2 (x z)/(x^2 + 2 x y + y^2) + z^2/(x^2 + 2 x y + y^2)

In[3]:= ExpandAll[E^(I a (t-b))]
Out[3]= E^(-I a b + I a t)

In[4]:= ExpandAll[((1+a) (1+b))[x]]
Out[4]= (1 + a + b + a b)[x]

In[5]:= ExpandAll[(f[(x+y)^2] + g[(y+z)^2])^2, x]
Out[5]= f[x^2 + 2 x y + y^2]^2 + g[(y + z)^2]^2 + 2 f[x^2 + 2 x y + y^2] g[(y + z)^2]
```

### Applications (4)

Reaches inside a function argument, unlike Expand

```mathematica
In[6]:= ExpandAll[Sin[(x + 1)^2]]
Out[6]= Sin[1 + 2 x + x^2]
```

The exponent is expanded too

```mathematica
In[7]:= ExpandAll[Exp[(a + b)^2]]
Out[7]= E^(a^2 + 2 a b + b^2)
```

The denominator is expanded as well

```mathematica
In[8]:= ExpandAll[(x + 1)^2/(x + 2)^2]
Out[8]= 1/(4 + 4 x + x^2) + 2 x/(4 + 4 x + x^2) + x^2/(4 + 4 x + x^2)
```

A second argument confines expansion to parts containing x

```mathematica
In[9]:= ExpandAll[(x + 1)^2 (y + 1)^2, x]
Out[9]= (1 + y)^2 + 2 x (1 + y)^2 + x^2 (1 + y)^2
```

## Implementation notes

**Algorithm.** `builtin_expand_all` drives the recursive `expr_expand_all_impl`.
Where `Expand` distributes only at the top level, `ExpandAll` reaches every
subexpression — function heads and arguments, exponents, and the bases of
denominators — expanding bottom-up: it recurses into the head and each argument
first (an `Inequality`'s operator-symbol slots at odd indices are passed through
untouched), rebuilds the node, and then applies the top-level distributor
`expr_expand_impl` at that level. A denominator factor `Power[base, -k]`
(`k > 0`) is handled specially, mirroring `ExpandDenominator`: recurse into the
`base`, expand it to the positive power, and keep the reciprocal, so a
surrounding `Times` later distributes the numerator across the expanded
denominator. `ExpandAll[expr, patt]` leaves any part free of `patt` unchanged
(`expr_contains_patt` guards each descent). Like `Expand`, an expansion too large
to fit in memory returns `Overflow[]` rather than silently declining
(`overflow_mode` is set on the user-facing entry).

**Data structures.** Pure `Expr`-tree recursion; each level allocates a fresh
`Expr**` argument array, rebuilds through `eval_and_free`, and hands the result
to `expr_expand_impl` (which itself uses the FLINT polynomial multiplier for
polynomial-over-`Q` factors). No persistent state beyond the recursion stack.

**Complexity / limits.** Bounded by the size of the fully-expanded tree; the
`Overflow[]` guard uses the same Newton-box size estimate as `Expand`. A symbolic
structural head — no packed/NDArray or `Compile[]` path. Threads over equations,
inequalities, logic functions, and lists; `Protected`.

- `Protected`.
- A thin recursive driver over the accelerated `Expand`: it descends into every
  part of `expr` — function heads, arguments, exponents, and the bases of
  denominators — and applies `Expand` (and `ExpandDenominator`) at each node,
  where a top-level `Expand` reaches none of them.
- `ExpandAll[expr, patt]` avoids expanding parts of `expr` that do not contain
  terms matching the pattern `patt`.
- Threads over `List`, equations, inequalities, and logic functions.
- Called with a number of arguments other than 1 or 2, emits `ExpandAll::argt`
  and stays unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [Expand](../../algebra/Expand/), [ExpandDenominator](../../algebra/ExpandDenominator/), [List](../../other-advanced/List/)

- Source: [`src/expand.c`](https://github.com/stblake/mathilda/blob/main/src/expand.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)
- Tests: [`tests/test_expand.c`](https://github.com/stblake/mathilda/blob/main/tests/test_expand.c)

## Notes & additional examples

### Notes

`ExpandAll[expr]` expands products and integer powers in *every* part of `expr` —
function heads and arguments, exponents, and the bases of denominators — where a
plain `Expand` distributes only at the top level. So `ExpandAll[Sin[(x + 1)^2]]`
expands the argument to `Sin[1 + 2 x + x^2]`, which `Expand` leaves untouched, and
`ExpandAll[(x + 1)^2/(x + 2)^2]` expands the `(x + 2)^2` denominator. The
two-argument form `ExpandAll[expr, patt]` leaves any part free of `patt` alone.
Like `Expand`, an expansion too large to fit in memory returns `Overflow[]` rather
than declining; it threads over lists, equations, inequalities, and logic
functions.
