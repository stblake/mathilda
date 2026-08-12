# RootReduce

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RootReduce[expr] canonicalises an algebraic expression: a constant algebraic number becomes a rational, a quadratic radical, or a Root object; a rational function over a radical tower has its denominator rationalised; a polynomial/rational function in a free variable has its constant-algebraic coefficients canonicalised. Threads over lists, equations, inequalities and logic. Option: Method -> "Automatic" | "Recursive" | "NumberField".`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (8)

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

In[6]:= RootReduce[1/(1 + k^(1/3))]        (* parametric tower *)
Out[6]= (1 - k^(1/3) + k^(2/3))/(1 + k)

In[7]:= RootReduce[(Sqrt[2] + Sqrt[3] - Sqrt[5 + 2 Sqrt[6]]) x^2 + x + 1]
Out[7]= 1 + x

In[8]:= RootReduce[a x^2 + Sqrt[8] x]      (* thread over coefficients *)
Out[8]= 2 Sqrt[2] x + a x^2
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

RootReduce also threads over equations, inequalities and logic functions (Equal, Less, And, ...), and for equations/inequalities of constant algebraic numbers it decides the (in)equality exactly via `qqbar`. It is Listable, so it threads over lists elementwise.

Options: Method -> "Automatic" | "Recursive" | "NumberField" (see flint_qqbar).

When `expr` carries no algebraic content (or the case is out of scope) it is returned unchanged, matching WL. Ownership follows the builtin contract: return a new tree or steal from `res`; never expr_free(res).

## Implementation notes

- `Protected`, `Listable`. Threads over lists, and over equations, inequalities
  and logic functions (`Equal`, `Unequal`, `Less`, `And`, ...); for
  (in)equalities of constant algebraic numbers it decides the relation exactly
  via `qqbar`.
- `Method`: `"Recursive"`/`"Automatic"` fold `qqbar` arithmetic bottom-up;
  `"NumberField"` re-expresses the value through a single primitive element of a
  common number field (`qqbar_express_in_field`). All three yield the identical
  canonical result. A `Root[]` object of degree ≤ 2 (or degree 1) auto-reduces
  to a quadratic radical / rational.
- One positional argument is required; other arg counts emit `RootReduce::argx`.
  An unknown `Method` emits `RootReduce::mtd`. Idempotent.

**Attributes:** `Listable`, `Protected`.

## See also

[Power](../../arithmetic/Power/), [Root](../../solutions-of-equations/Root/), [Re](../../arithmetic/Re/), [Im](../../arithmetic/Im/), [Cancel](../../algebra/Cancel/), [Equal](../../comparisons/Equal/), [Unequal](../../comparisons/Unequal/), [Less](../../comparisons/Less/)

## References

- Source: [`src/rootreduce.c`](https://github.com/stblake/mathilda/blob/main/src/rootreduce.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
- Tests: [`tests/test_rootreduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_rootreduce.c)
