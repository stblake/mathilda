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

## Examples (13)

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
Out[11]= {{x -> 0.000903027}}

In[12]:= FindInstance[a^2 + b c == 0 && a b + b d == 0 && a c + c d == 0 && b c + d^2 == 0 && a d - b c != 0, {a, b, c, d}, Reals]
Out[12]= {}

In[13]:= FindInstance[Xor[p, q] && Implies[q, r] && Not[Equivalent[p, r]], {p, q, r}, Booleans]
Out[13]= {{p -> True, q -> False, r -> False}}
```

## Algorithm

reduce_companions.c

```text
Companion builtins for `Reduce` (REDUCE_PLAN.md, Phase 8).  v1: LogicalExpand
+ a minimal NotElement head.
```

LogicalExpand distributes a logical statement to disjunctive normal form (an Or of Ands of literals), applying idempotence / complementation / absorption contractions, and collapsing to True (tautology) or False (contradiction)

```text
when the statement decides.  Every non-connective subexpression is treated as
```

an OPAQUE Boolean atom -- a symbol, a relation `x == a`, a membership

```text
`Element[..]` -- with NO domain reasoning, exactly as Mathematica's
LogicalExpand does.  Two relational atoms are complementary iff one is the
```

(head-flipped) logical negation of the other (`x==a` / `x!=a`, `x<1` / `x>=1`,

```text
`Element` / `NotElement`, `a` / `!a`).
```

The True/False collapse is sound *and* complete without truth-table enumeration: over independent opaque atoms a DNF is unsatisfiable iff every clause holds a complementary pair -- i.e. it distributes to ZERO surviving

```text
clauses.  So `phi` empty => False, and `Not[phi]` empty => True (the negation
```

is unsatisfiable, hence `phi` is a tautology).

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Reduce](../../solutions-of-equations/Reduce/), [Solve](../../solutions-of-equations/Solve/), [LogicalExpand](../../solutions-of-equations/LogicalExpand/), [E](../../mathematical-constants/E/), [Tan](../../elementary-functions/Tan/), [Sec](../../elementary-functions/Sec/), [Cot](../../elementary-functions/Cot/), [Csc](../../elementary-functions/Csc/)

- Source: [`src/solve/reduce_companions.c`](https://github.com/stblake/mathilda/blob/main/src/solve/reduce_companions.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)
