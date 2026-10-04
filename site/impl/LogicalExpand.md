---
source: src/solve/reduce_companions.c
---
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
