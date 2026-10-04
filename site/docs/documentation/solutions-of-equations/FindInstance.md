# FindInstance

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindInstance[expr, vars]`**

Finds a single instance of vars satisfying the statement expr -- a logical combination of equations and inequalities -- returned in Solve's form {{x -\> v, ...}}, or {} if none exists.  The default domain is Complexes, or Reals when expr carries an ordering (as in Reduce).

**`FindInstance[expr, vars, dom]`**

Finds an instance over dom: Complexes, Reals, Integers, Rationals, or Booleans (Boolean satisfiability).

**`FindInstance[expr, vars, dom, n]`**

Finds up to n instances (fewer if fewer exist).

<details>
<summary>Notes</summary>

Every instance returned is verified against expr, so it is always a true solution.  Variables may be symbols or indexed forms c\[i\]. FindInstance may find an instance even where Reduce cannot give a complete reduction -- instantiating parametric Diophantine families, searching a bounded integer box over the Integers, and (for transcendental or inexact Real systems) a numerical feasibility search. It returns {} only when the set is provably empty -- including a Groebner certificate for declined polynomial systems -- and stays unevaluated otherwise.  Modulus -\> p over Z/pZ.

</details>

## Examples (19)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (13)

```mathematica
In[1]:= FindInstance[x^2 == 2, x]
Out[1]= {{x -> -Sqrt[2]}}

In[2]:= FindInstance[x^2 + y^2 <= 1, {x, y}, Reals]
Out[2]= {{x -> 0, y -> 0}}

In[3]:= FindInstance[x^2 - 3 y^2 == 1 && 10 < x < 100, {x, y}, Integers]
Out[3]= {{x -> 26, y -> -15}}

In[4]:= FindInstance[x^2 < 10 && x > 0, x, Integers, 3]
Out[4]= {{x -> 1}, {x -> 2}, {x -> 3}}

In[5]:= FindInstance[Xor[a, b, c, d] && (a || b) && ! (c || d), {a, b, c, d}, Booleans]
Out[5]= {{a -> True, b -> False, c -> False, d -> False}}

In[6]:= FindInstance[x^2 + y^3 == 3 && x + 2 y >= 4 && x y == 5, {x, y}, Reals]
Out[6]= {}

In[7]:= FindInstance[x^2 - 61 y^2 == 1 && x > 0 && y > 0, {x, y}, Integers]
Out[7]= {{x -> 1766319049, y -> 226153980}}

In[8]:= FindInstance[Sin[1/x] == 0 && 0 < x < 10^-5, x, Reals]
Out[8]= {{x -> 1/31831/Pi}}

In[9]:= FindInstance[a^3 + b^3 + c^3 == d^3 && a > 0 && b > 0 && c > 0 && d > 0, {a, b, c, d}, Integers]
Out[9]= {{a -> 5, b -> 4, c -> 3, d -> 6}}

In[10]:= FindInstance[Total[Array[c, 15]*Prime[Range[15]]] == 500 && And @@ Thread[0 <= Array[c, 15] <= 1], Array[c, 15], Integers]
Out[10]= {}

In[11]:= FindInstance[0 < x < 0.001 && Sin[1/x] > 0.999, x, Reals]
Out[11]= {{x -> 0.000827826}}

In[12]:= FindInstance[a^2 + b c == 0 && a b + b d == 0 && a c + c d == 0 && b c + d^2 == 0 && a d - b c != 0, {a, b, c, d}, Reals]
Out[12]= {}

In[13]:= FindInstance[Xor[p, q] && Implies[q, r] && Not[Equivalent[p, r]], {p, q, r}, Booleans]
Out[13]= {{p -> True, q -> False, r -> False}}
```

### Applications (6)

One witness to a linear system

```mathematica
In[14]:= FindInstance[x + y == 10 && x - y == 2, {x, y}]
Out[14]= {{x -> 6, y -> 4}}
```

A single root of a quadratic over the complexes

```mathematica
In[15]:= FindInstance[x^2 == 2, x]
Out[15]= {{x -> -Sqrt[2]}}
```

A point on the unit circle in the open first quadrant

```mathematica
In[16]:= FindInstance[x^2 + y^2 == 1 && x > 0 && y > 0, {x, y}, Reals]
Out[16]= {{x -> 1/2, y -> 1/2 Sqrt[3]}}
```

An integer lattice point on the line

```mathematica
In[17]:= FindInstance[2 x + 3 y == 1, {x, y}, Integers]
Out[17]= {{x -> -1, y -> 1}}
```

Provably empty: returns {}

```mathematica
In[18]:= FindInstance[x^2 + 1 == 0 && x > 0, x, Reals]
Out[18]= {}
```

A satisfying Boolean assignment

```mathematica
In[19]:= FindInstance[a && Xor[a, b], {a, b}, Booleans]
Out[19]= {{a -> True, b -> False}}
```

## Implementation notes

**Algorithm.** `builtin_find_instance` peels trailing options (`Modulus` is
honoured; `Method` / `WorkingPrecision` / `RandomSeeding` are accepted and
ignored, since the search is exact and deterministic), reads the optional domain
and witness count `n`, then runs `fi_run_search` under a message mute — its
internal `Reduce` / `Solve` / `NMinimize` probes are speculative and, as in
Mathematica, must not leak diagnostics. The strategy is **soundness-first and
verify-gated**: every candidate point is checked against the original statement
(`expr /. point === True`, by `fi_verify`) before being accepted, so a returned
instance is never wrong, and a point that cannot be verified is discarded rather
than reported.

`fi_run_search` is a cascade of witness sources, tried cheapest/most-exact first.
`Reduce` is the satisfiability-and-solution-set oracle (step 1): its `False` is
taken as `{}` only for exactly-decidable systems — for transcendental or inexact
systems (`fi_is_transc_inexact`) that `False` is distrusted, and indexed variables
`c[i]` skip the oracle since `Reduce` rejects them — otherwise each top-level `Or`
clause is sampled into a point. Then `Solve` fallbacks with generated-parameter
instantiation reach parametric Diophantine families (the Pell equation); a bounded
integer box search runs over `Integers`; an equations-only retry, a 1-variable real
transcendental root bracket, a solve-one-then-sample grid, and Rabinowitsch ideal
saturation each cover systems the oracles decline; a Gröbner emptiness certificate
proves `{}` for a declined polynomial system; structured exact candidate sampling
finds branch-cut and open-region witnesses; and numerical feasibility
(`NMinimize` / least-infeasibility) is the last resort for transcendental/inexact
`Real` systems. Interval samples come from `rru_rational_between`. The `Booleans`
domain is handled separately by `fi_boolean`, reusing the `LogicalExpand` DNF
engine for satisfiability.

**Data structures.** Witnesses accumulate in a `FiWit` buffer and are returned in
`Solve`'s rule-list form `{{x -> v, ...}, ...}`. A "variable" is matched
structurally with `expr_eq`, so a plain symbol (`x`) and an indexed form (`c[i]`)
are both accepted. Everything is `Expr` trees driven through `evaluate`.

**Complexity / limits.** The default domain is `Complexes`, or `Reals` when the
statement carries an ordering (as in `Reduce`). `{}` is returned only when the set
is *provably* empty — a `Reduce` `False` on an exactly-decidable system, an
exhausted finite integer box, or a Gröbner certificate; when no witness is found
and emptiness is not proved, the call stays unevaluated. Cost is dominated by the
`Reduce`/`Solve`/Gröbner probes the cascade invokes; `Modulus -> p` searches over
`Z/pZ`.

**Attributes:** `Protected`.

## References

**See also:** [Reduce](../../solutions-of-equations/Reduce/), [Solve](../../solutions-of-equations/Solve/), [LogicalExpand](../../solutions-of-equations/LogicalExpand/), [E](../../mathematical-constants/E/), [Tan](../../elementary-functions/Tan/), [Sec](../../elementary-functions/Sec/), [Cot](../../elementary-functions/Cot/), [Csc](../../elementary-functions/Csc/)

- Source: [`src/solve/reduce_companions.c`](https://github.com/stblake/mathilda/blob/main/src/solve/reduce_companions.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)

## Notes & additional examples

### Notes

`FindInstance[expr, vars]` returns *one* instance of `vars` that satisfies the
statement `expr`, in `Solve`'s rule-list form `{{x -> v, ...}}`; `{}` means the
solution set is provably empty. `FindInstance[expr, vars, dom]` names the domain
(`Complexes`, `Reals`, `Integers`, `Rationals`, or `Booleans`), and a trailing
integer `n` asks for up to `n` instances. The default domain is `Complexes`, or
`Reals` when `expr` carries an ordering — the same rule `Reduce` uses.

**Every instance returned is verified** against the original statement, so a
reported point is always a true solution. Because of that, `FindInstance` can
succeed where `Reduce` gives no complete reduction: it instantiates parametric
Diophantine families, searches a bounded integer box over `Integers`, finds
branch-cut and open-region witnesses by structured sampling, and falls back to a
numerical feasibility search for transcendental or inexact real systems.
Variables may be plain symbols or indexed forms `c[i]`.

A returned `{}` is a genuine emptiness proof — from a decidable `Reduce`, an
exhausted finite integer box, or a Gröbner certificate — whereas a system where
no witness is found and emptiness is not proved stays unevaluated. Use
`Modulus -> p` to search over `Z/pZ`.
