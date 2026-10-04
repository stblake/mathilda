# RootReduce

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RootReduce[expr] canonicalises an algebraic expression: a constant algebraic number becomes a rational, a quadratic radical, or a Root object; a rational function over a radical tower has its denominator rationalised; a polynomial/rational function in a free variable has its constant-algebraic coefficients canonicalised. Threads over lists, rules (Solve results), equations, inequalities and logic. Option: Method -> "Automatic" | "Recursive" | "NumberField".`**

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (9)

```mathematica
In[1]:= RootReduce[Sqrt[2] + Sqrt[3]]
Out[1]= Root[1 - 10 #1^2 + #1^4 &, 4]

In[2]:= RootReduce[(Sqrt[18] + Sqrt[27]) / Sqrt[5 + 2 Sqrt[6]]]
Out[2]= 3

In[3]:= RootReduce[1/(1 + Sqrt[2])]
Out[3]= -1 + Sqrt[2]

In[4]:= RootReduce[1/(1 + 2^(1/3) + 2^(2/3))]
Out[4]= Root[-1 + 3 #1 + 3 #1^2 + #1^3 &, 1]

In[5]:= RootReduce[Sqrt[2] + Sqrt[3] + Sqrt[5] == Sqrt[10 + 2 Sqrt[15] + 4 Sqrt[4 + Sqrt[15]]]]
Out[5]= True
```

Parametric tower

```mathematica
In[6]:= RootReduce[1/(1 + k^(1/3))]
Out[6]= (1 - k^(1/3) + k^(2/3))/(1 + k)
```

```mathematica
In[7]:= RootReduce[(Sqrt[2] + Sqrt[3] - Sqrt[5 + 2 Sqrt[6]]) x^2 + x + 1]
Out[7]= 1 + x
```

Thread over coefficients

```mathematica
In[8]:= RootReduce[a x^2 + Sqrt[8] x]
Out[8]= 2 Sqrt[2] x + a x^2
```

Thread over Solve rules

```mathematica
In[9]:= RootReduce[{u -> Sqrt[8], v -> 1/(1 + Sqrt[2])}]
Out[9]= {u -> 2 Sqrt[2], v -> -1 + Sqrt[2]}
```

### Applications (6)

A degree-4 algebraic number: a Root object

```mathematica
In[10]:= RootReduce[Sqrt[2] + Sqrt[3]]
Out[10]= Root[1 - 10 #1^2 + #1^4 &, 4]
```

Pulled apart into 2 Sqrt[2]

```mathematica
In[11]:= RootReduce[Sqrt[8]]
Out[11]= 2 Sqrt[2]
```

Combined into Sqrt[6]

```mathematica
In[12]:= RootReduce[Sqrt[2] Sqrt[3]]
Out[12]= Sqrt[6]
```

A nested radical, denested

```mathematica
In[13]:= RootReduce[Sqrt[3 + 2 Sqrt[2]]]
Out[13]= 1 + Sqrt[2]
```

Decided to be exactly zero

```mathematica
In[14]:= RootReduce[(1 + Sqrt[5])/2 - GoldenRatio]
Out[14]= 0
```

An equality, decided exactly

```mathematica
In[15]:= RootReduce[Sqrt[2] + Sqrt[3] == Sqrt[5 + 2 Sqrt[6]]]
Out[15]= True
```

## Algorithm

Mathilda — RootReduce implementation.

RootReduce[expr] canonicalises an algebraic expression. It dispatches between two rigorous FLINT engines depending on the shape of `expr`:

```text
  (1) Constant algebraic NUMBERS (no free symbol) — integers, rationals,
      radicals, roots of unity, the imaginary unit and Root[] objects
      combined by +,-,*,/,^ — are canonicalised via FLINT `qqbar`
      (src/poly/flint_qqbar.c) to a single representative: a rational, a
      quadratic radical expression, or a Root[Function[minpoly&], k] object.
      This is WL's central RootReduce behaviour.

  (2) Algebraic FUNCTIONS over a tower Q(params)(radicals) — radicals whose
      radicand carries a free variable (e.g. the Goursat k^(1/3) towers) —
      are rationalised by flint_algebraic_field_canonical (src/poly/
      flint_bridge.c): the denominator is inverted in the field by an exact
      linear solve, no numeric oracle.

  (3) POLYNOMIALS / RATIONAL FUNCTIONS in a free variable whose coefficients
      are constant algebraic numbers — threaded over via rr_thread_coeffs:
      each maximal constant-algebraic subexpression (a coefficient) is
      canonicalised via qqbar and the free-variable structure is left intact,
      so a vanishing radical coefficient reduces to 0 and its monomial drops
      out. Plain polynomial cancellation is NOT performed (that is Cancel).
```

RootReduce also threads over equations, inequalities and logic functions (Equal, Less, And, ...), and for equations/inequalities of constant algebraic numbers it decides the (in)equality exactly via `qqbar`. It threads over an (immediate) Rule too, so a Solve result {u -> value, ...} reduces the same way the corresponding Reduce result does. It is Listable, so it threads over lists elementwise.

Options: Method -> "Automatic" | "Recursive" | "NumberField" (see flint_qqbar).

When `expr` carries no algebraic content (or the case is out of scope) it is returned unchanged, matching WL. Ownership follows the builtin contract: return a new tree or steal from `res`; never expr_free(res).

## Implementation notes

**Algorithm.** `builtin_rootreduce` is a dispatcher over three rigorous FLINT
engines, chosen by the shape of the argument. It first splits a trailing
`Method -> "Automatic" | "Recursive" | "NumberField"` option (recognised by
reading `Options[RootReduce]`, so a `Solve`-result rule `u -> value` is kept as
a positional argument, not mistaken for an option).

1. **Constant algebraic number** (no free symbol) → `flint_qqbar_canonical`
   (`src/poly/flint_qqbar.c`). `to_qqbar` converts an expression built from
   integers, rationals, radicals `Power[base, p/q]`, roots of unity `(−1)^(p/q)`,
   the imaginary unit, `Root[]` objects and `GoldenRatio`, combined by
   `+ − * / ^`, into a single `qqbar_t`. `qqbar_to_expr` then renders a canonical
   representative: a **rational** when `qqbar_is_rational`; `(a + b Sqrt[c])/q`
   for **degree 2** via `qqbar_get_quadratic`; and `Root[Function[minpoly &], k]`
   for **degree ≥ 3**, with the index `k` placed by WL's root ordering
   (`wl_root_index`) and the whole object memoised by `(minpoly, k)`.
   `Method -> "NumberField"` re-expresses the value through a single primitive
   element of the field (`number_field_value`).
2. **Parametric algebraic function** — a radical whose radicand carries a free
   variable — → `flint_algebraic_field_canonical` (`src/poly/flint_bridge.c`):
   the denominator is inverted in the field by an exact linear solve over the
   tower, with no numeric oracle.
3. **Polynomial / rational function in a free variable** with
   constant-algebraic coefficients → `flint_qqbar_reduce_coeffs`: each maximal
   constant-algebraic subexpression (a coefficient) is canonicalised via `qqbar`
   while the free-variable structure is left intact, so a vanishing radical
   coefficient reduces to `0` and its monomial drops out. (Plain polynomial
   cancellation is `Cancel`'s job, not done here.)

`RootReduce` also **threads** over `Equal`/`Unequal`/`Less`/`…`/`And`/`Or` and
over an immediate `Rule` (a `Solve` result entry). For a binary (in)equality of
constant algebraic numbers it is *decided exactly* — `flint_qqbar_equal`
(1/0/−1) for `Equal`/`Unequal`, `flint_qqbar_compare` (sign, or −2 undecided)
for the ordering relations — so `Sqrt[2] + Sqrt[3] == Sqrt[5 + 2 Sqrt[6]]`
returns `True` with no floating point. When the argument carries no algebraic
content it is returned unchanged (the positional arg is stolen out of `res`).

**Data structures.** FLINT `qqbar_t` (an exact algebraic number: primitive
integer minimal polynomial as an `fmpz_poly` plus an isolating complex
enclosure). Three session caches accelerate the hot paths: a `to_qqbar`
conversion cache, the shared `wl-roots` ordering cache, and a `(minpoly, k)` →
`Root[]`-object memo (`q2e_cache`). Equality and comparison are rigorous
decisions on the minimal polynomials, **not** a numeric zero test.

**Complexity / limits.** Dominated by `qqbar` arithmetic over the compositum of
the atoms; bounded by the degree cap `QQBAR_DEGREE_CAP = 120`, above which it
declines. A non-real ordering comparison is undecided (`−2`) and leaves the
relation unevaluated rather than guessing. `Listable`, `Protected`; option
`Method -> "Automatic" | "Recursive" | "NumberField"`.

- `Protected`, `Listable`. Threads over lists, over equations, inequalities and
  logic functions (`Equal`, `Unequal`, `Less`, `And`, ...), and over an
  (immediate) `Rule` — so `Solve[...] // RootReduce` reduces the right-hand side
  of each `u -> value` entry the same way `Reduce[...] // RootReduce` reduces
  each `u == value`, leaving the free-variable left-hand side intact. For
  (in)equalities of constant algebraic numbers it decides the relation exactly
  via `qqbar`. A `Rule` whose left-hand side is the option name `Method` is a
  trailing option; any other symbol left-hand side (e.g. `u -> value`) is a
  positional argument threaded over, not an option.
- `Method`: `"Recursive"`/`"Automatic"` fold `qqbar` arithmetic bottom-up;
  `"NumberField"` re-expresses the value through a single primitive element of a
  common number field (`qqbar_express_in_field`). All three yield the identical
  canonical result. A `Root[]` object of degree ≤ 2 (or degree 1) auto-reduces
  to a quadratic radical / rational.
- One positional argument is required; other arg counts emit `RootReduce::argx`.
  An unknown `Method` emits `RootReduce::mtd`. Idempotent.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [Power](../../arithmetic/Power/), [Root](../../solutions-of-equations/Root/), [Re](../../arithmetic/Re/), [Im](../../arithmetic/Im/), [Cancel](../../algebra/Cancel/), [Equal](../../comparisons/Equal/), [Unequal](../../comparisons/Unequal/), [Less](../../comparisons/Less/)

- H. Cohen, *A Course in Computational Algebraic Number Theory*, GTM 138 (Springer, 1993), ch. 4 (algebraic numbers, minimal polynomials, number fields).
- The FLINT library (https://flintlib.org), `qqbar` module — exact real and complex algebraic numbers via minimal polynomial plus isolating enclosure, with no numeric zero oracle.
- Source: [`src/rootreduce.c`](https://github.com/stblake/mathilda/blob/main/src/rootreduce.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
- Tests: [`tests/test_algebraicnumber.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumber.c)
- Tests: [`tests/test_algebraicnumberpolynomial.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumberpolynomial.c)
- Tests: [`tests/test_nf_rowreduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nf_rowreduce.c)
- Tests: [`tests/test_nullspace.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nullspace.c)

## Notes & additional examples

### Notes

`RootReduce[expr]` canonicalises an algebraic expression to a single
representative. A constant algebraic number becomes a rational, a quadratic
radical `(a + b Sqrt[c])/q`, or a `Root[poly &, k]` object for degree three and
up; the representative is unique, so two spellings of the same number reduce to
the same thing and their difference reduces to `0`.

Reduction runs on FLINT's exact `qqbar` engine — minimal polynomial plus an
isolating enclosure, with **no numeric zero oracle** — so `RootReduce` also
*decides* equations and inequalities between constant algebraic numbers exactly,
and threads over lists, `Solve`-result rules, and logical combinations. It
leaves anything with no algebraic content unchanged. Option
`Method -> "Automatic" | "Recursive" | "NumberField"`.
