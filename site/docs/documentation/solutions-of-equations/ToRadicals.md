# ToRadicals

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ToRadicals[expr]`**

attempts to express all Root objects in expr in terms of radicals.

<details>
<summary>Notes</summary>

ToRadicals can always give expressions in terms of radicals when the highest degree of the polynomial that appears in any Root object is four.  Binomial Root objects of the form Root\[Function\[a #^n + b\], k\] are also reduced to radicals for any degree n.  Other Root objects of degree five or higher are returned unchanged. If Root objects in expr contain parameters, ToRadicals\[expr\] may yield a result that is not equal to expr for all values of the parameters. ToRadicals automatically threads over lists, equations, inequalities, and logic functions.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= ToRadicals[Root[Function[#^2 + 3 # + 5], 1]]
Out[1]= 1/2 (-3 - I Sqrt[11])

In[2]:= ToRadicals[Root[Function[#^2 + 3 # + 5], 2]]
Out[2]= 1/2 (-3 + I Sqrt[11])

In[3]:= ToRadicals[Root[Function[#^5 - 2], 3]]
Out[3]= (-1)^(4/5) 2^(1/5)

In[4]:= With[{r = ToRadicals[Root[Function[#^4 + 3 #^3 - 5 #^2 - 7 # + 9], 1]]}, Chop[N[r^4 + 3 r^3 - 5 r^2 - 7 r + 9, 30]]]
Out[4]= 0

In[5]:= ToRadicals[Root[Function[#^5 - # - 1], 1]]      (* non-binomial deg 5 *)
Out[5]= Root[#1^5 - #1 - 1 &, 1]

In[6]:= ToRadicals[Root[Function[#^2 - 2], 2] < 3]      (* threading *)
Out[6]= True
```

### Applications (4)

```mathematica
In[1]:= ToRadicals[Root[#^2 - 2 &, 1]]
Out[1]= -Sqrt[2]
```

The smaller root of `x^2 + x - 1` is the negative reciprocal of the golden ratio,
recovered exactly in radicals:

```mathematica
In[1]:= ToRadicals[Root[#^2 + # - 1 &, 1]]
Out[1]= 1/2 (-1 - Sqrt[5])
```

Cardano's formula appears automatically for the real root of a depressed cubic:

```mathematica
In[1]:= ToRadicals[Root[1 + #1 + #1^3 &, 1]]
Out[1]= -1/3 ((1/2 (27 + 3 Sqrt[93]))^(1/3) - 3/(1/2 (27 + 3 Sqrt[93]))^(1/3))
```

It threads through the implicit `Root` objects produced by `Solve`, here giving
the three cube roots of two:

```mathematica
In[1]:= ToRadicals[Solve[x^3 - 2 == 0, x]]
Out[1]= {{x -> 2^(1/3)}, {x -> -(-1)^(1/3) 2^(1/3)}, {x -> (-1)^(2/3) 2^(1/3)}}
```

## Algorithm

radicals.c

ToRadicals: convert held Root[Function[poly], k] objects into closed-form

```text
radical expressions.  See radicals.h for the public contract.
```

Algorithm (per Root node):

```text
  1. Extract the polynomial body from Root[Function[..], k].  Both the
     Slot[1] form `Function[expr]` and the bound-variable form
     `Function[t, expr]` are accepted.
  2. Substitute Slot[1] (or t) with a fresh symbol `x$` so the existing
     get_coeff / get_degree_poly polynomial machinery operates on a
     standard univariate polynomial.
  3. Dispatch on degree d:
       d == 1 : linear, x = -c0/c1
       d == 2 : quadratic formula
       d == 3 : Cardano
       d == 4 : Ferrari (depressed quartic + resolvent cubic)
       d >= 5 : binomial fast-path a x^n + b only; otherwise leave the
                Root untouched.
     Each path produces ALL d radical roots as a freshly-owned Expr**.
  4. Select the k-th root in Mathilda's canonical Root ordering by
     computing N[Root[poly, k]] at machine precision (via
     root_numericalize) and picking the radical root whose numeric
     value lies closest in the complex plane.
     When numeric evaluation is unavailable (the polynomial has
     parametric coefficients), fall back to the natural per-formula
     order with k - 1 as the index.
```

Threading: the top-level walker is a structural recurrence that reconstructs every EXPR_FUNCTION node it visits, so a Root buried inside List, Equal, Less, And, Or, ... is processed identically.

Memory: every internal helper returns a freshly-owned Expr*; inputs

```text
are borrowed and deep-copied wherever they appear in the output.  The
```

exact-arithmetic core (eval_and_free, Plus/Times/Power normalisation) does the bookkeeping for intermediate trees.

## Implementation notes

**Algorithm.** `builtin_to_radicals` (`builtin_to_radicals`) converts held `Root[Function[poly], k]` objects into closed-form radical expressions. The top-level walker is a structural recurrence that rebuilds every `EXPR_FUNCTION` node, so a `Root` buried inside `List`/`Equal`/`Less`/`And`/`Or`/... is handled identically.

Per `Root` node: (1) extract the polynomial body, accepting both the `Slot[1]` form `Function[expr]` and the bound-variable form `Function[t, expr]`; (2) substitute the slot/variable with a fresh symbol `x$` so the standard `get_coeff`/`get_degree_poly` univariate machinery applies; (3) dispatch on degree d — `d=1` linear (`-c0/c1`), `d=2` quadratic formula, `d=3` Cardano, `d=4` Ferrari (depressed quartic + resolvent cubic), `d≥5` only the binomial fast-path `a·x^n + b`, otherwise the `Root` is left untouched — each path producing all d radical roots as a fresh `Expr**`; (4) select the k-th root in Mathilda's canonical `Root` ordering by computing `N[Root[poly, k]]` at machine precision (`root_numericalize`) and picking the radical root closest in the complex plane, falling back to the natural per-formula order (index `k-1`) when coefficients are parametric and numeric evaluation is unavailable.

**Data structures.** `Expr*` trees; degree dispatch reuses the polynomial coefficient extractors from `src/poly/poly.c`. Intermediate radical-expression bookkeeping rides on `eval_and_free` and the Plus/Times/Power normalisation. Inputs are borrowed and deep-copied into the output.

**Complexity / limits.** Closed forms exist only up to degree 4 (Abel–Ruffini); degree ≥ 5 is supported solely for binomials. Root selection costs one numeric `Root` evaluation per node.

- `Protected`.
- Closed-form radicals are always returned when the polynomial has degree
  at most four — linear (trivial), quadratic (`Sqrt`), cubic (Cardano), and
  quartic (Ferrari via the depressed quartic + resolvent cubic).
- Binomial Root objects `Root[Function[a #^n + b], k]` are reduced to
  radicals for any degree `n`, using the principal `n`-th root multiplied
  by `(-1)^(2 (k-1) / n)`.
- Other Root objects of degree ≥ 5 are returned unchanged — the system
  makes no attempt at decomposition or solvable-Galois detection (cf.
  Mathematica's note "ToRadicals cannot find them").
- The k-th radical root is selected to agree with `N[Root[poly, k]]`'s
  canonical ordering (real-first ascending, complex by `Re` / `|Im|` /
  negative-`Im` first) — each formula's natural emission order is
  numerically matched against `root_numericalize` at machine precision.
  When the polynomial carries parametric coefficients (no numericalisation
  possible), the natural per-formula index `k - 1` is used and the result
  is allowed to disagree with `expr` for some parameter values, matching
  Mathematica's `nongen` behaviour.
- Walks its argument recursively, so `Root[..]` nodes inside `List`,
  `Equal`, `Less`, `Greater`, `And`, `Or`, `Not`, `Implies`, ... thread
  automatically — every `Root` anywhere in the tree is processed
  independently and the surrounding structure is preserved.
- Idempotent: `ToRadicals[ToRadicals[expr]] === ToRadicals[expr]`, since a
  successful conversion produces an expression free of `Root[..]` nodes.

**Attributes:** `Protected`.

## See also

[Sqrt](../../arithmetic/Sqrt/), [Re](../../arithmetic/Re/), [Im](../../arithmetic/Im/), [List](../../other-advanced/List/), [Equal](../../comparisons/Equal/), [Less](../../comparisons/Less/), [Greater](../../comparisons/Greater/), [Root](../../solutions-of-equations/Root/)

## References

- G. Cardano, *Ars Magna*, 1545 (cubic); L. Ferrari (quartic resolvent, via Cardano).
- Source: [`src/radicals.c`](https://github.com/stblake/mathilda/blob/main/src/radicals.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)
- Tests: [`tests/test_radicals.c`](https://github.com/stblake/mathilda/blob/main/tests/test_radicals.c)

## Notes & additional examples

### Notes

`ToRadicals[expr]` rewrites `Root` objects in `expr` using radicals. It always
succeeds when the underlying polynomial has degree at most four (and for
binomial `Root[a #^n + b &, k]` of any degree); degree-five-and-higher Root
objects are returned unchanged. It threads automatically over lists, equations,
and the results of `Solve`.
