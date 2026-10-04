# LogicalExpand

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LogicalExpand[expr]`**

Expands the logical combination expr -- of equations, inequalities and Boolean atoms -- into disjunctive normal form (an Or of Ands), applying distributive, De Morgan, idempotence, complementation and absorption laws, and expanding Implies and Xor.  Returns True for a tautology and False for a contradiction.  Every non-logical subexpression is treated as an opaque Boolean atom (no domain reasoning), so e.g. x==a and x!=a are complementary literals.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (6)

Distribute And over Or into disjunctive normal form

```mathematica
In[1]:= LogicalExpand[(a || b) && c]
Out[1]= a && c || b && c
```

Full distribution gives four conjunctive clauses

```mathematica
In[2]:= LogicalExpand[(a || b) && (c || d)]
Out[2]= a && c || a && d || b && c || b && d
```

Material implication becomes !a || b

```mathematica
In[3]:= LogicalExpand[Implies[a, b]]
Out[3]= Not[a] || b
```

Exclusive-or expands to its two odd-parity clauses

```mathematica
In[4]:= LogicalExpand[Xor[a, b]]
Out[4]= a && Not[b] || Not[a] && b
```

Equivalence is both-true-or-both-false

```mathematica
In[5]:= LogicalExpand[Equivalent[a, b]]
Out[5]= Not[a] && Not[b] || b && a
```

Complementation collapses to False

```mathematica
In[6]:= LogicalExpand[a && !a]
Out[6]= False
```

## Implementation notes

**Algorithm.** `builtin_logical_expand` takes a single logical statement and
rewrites it into **disjunctive normal form** — an `Or` of `And`s — by the mutually
recursive `to_dnf` / `to_dnf_neg`. `And` distributes across its arguments
(`dnf_and`), `Or` takes their union (`dnf_or`), `Not` flips to `to_dnf_neg` and
applies De Morgan, and `Implies`, `Xor` and `Equivalent` are expanded into their
`And`/`Or` definitions; an `Element` / `NotElement` over a container (a list or
`Alternatives`) distributes elementwise (`element_multi_dnf`). Clause construction
applies idempotence, complementation and absorption as it goes, so redundant and
contradictory literals are dropped.

The result is then classified. An empty DNF (`phi.n == 0`) is a contradiction, so
`False` is returned; a DNF carrying an empty clause is structurally satisfiable
everywhere, so `True`. Otherwise `LogicalExpand` also expands the negation and, if
`!e` is unsatisfiable (`nphi.n == 0`), reports the tautology `True`; failing that
it emits the `Or`-of-`And`s via `dnf_to_expr` and runs `evaluate` once to
canonicalise the `And`/`Or` ordering. Every non-logical subexpression is treated as
an **opaque Boolean atom** (`leaf_dnf`) with no domain reasoning — but `x == a` and
`x != a` (and `Element`/`NotElement` of the same arguments) are recognised as
complementary literals, so they cancel.

**Data structures.** The working form is a `Dnf`: an array of `Clause`s, each an
array of literal `Expr*`s, where a literal carries a negation flag so `a` and `!a`
are the same atom with opposite sign. `dnf_true` / `dnf_false` are the identity
elements of the two combinators, and `dnf_has_empty` is the tautology test.

**Complexity / limits.** DNF conversion is worst-case exponential in the number of
distinct atoms (the `And`-over-`Or` distribution is a cartesian product of
clauses). The procedure is purely structural Boolean algebra: it performs **no**
numeric or domain reasoning beyond the literal-complement recognition above, so the
truth of an atom like `x > 0` is never decided — that is `Reduce`'s job.

**Attributes:** `Protected`.

## References

**See also:** [Implies](../../control-flow/Implies/), [Xor](../../control-flow/Xor/), [Reduce](../../solutions-of-equations/Reduce/), [Element](../../simplification/Element/), [NotElement](../../solutions-of-equations/NotElement/)

- Source: [`src/solve/reduce_companions.c`](https://github.com/stblake/mathilda/blob/main/src/solve/reduce_companions.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)

## Notes & additional examples

### Notes

`LogicalExpand[expr]` rewrites a logical combination into **disjunctive normal
form** — an `Or` of `And`s. It distributes `And` over `Or`, applies De Morgan to
`Not`, and expands `Implies`, `Xor` and `Equivalent` into their `And`/`Or`
definitions, simplifying by idempotence, complementation and absorption as it
builds each clause.

A tautology collapses to `True` and a contradiction to `False`. Every non-logical
subexpression is treated as an opaque Boolean atom — there is **no domain
reasoning**, so `LogicalExpand` never decides whether `x > 0` is true — but a pair
like `x == a` and `x != a` (or `Element` and `NotElement` of the same arguments)
is recognised as complementary literals and cancels. For the full solution set of
equations and inequalities, with domain reasoning, use `Reduce`; `LogicalExpand`
is the purely structural Boolean normaliser.
