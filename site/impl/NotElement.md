---
source: src/solve/reduce_companions.c
---
**Algorithm.** `builtin_not_element` is the Boolean negation of `Element`. Given
`NotElement[x, dom]` it builds `Element[x, dom]`, evaluates it, and inverts the
verdict: an `Element` that decides to `True` yields `False`, one that decides to
`False` yields `True`. If `Element` does not resolve to a literal (it returns a
symbolic `Element[...]`, or anything that is not `True`/`False`), the builtin
returns `NULL`, so `NotElement[x, dom]` stays symbolic. All the membership
logic — the recognised domains (`Integers`, `Rationals`, `Reals`, `Complexes`,
`Algebraics`, `Booleans`, ...), and the distribution of a container first argument
(a `List` or `Alternatives`) over its members — lives in `Element`; `NotElement`
only flips the result.

In a logical statement fed to `Reduce`, `LogicalExpand`, or the quantifier engine,
`NotElement` is also understood directly by the DNF translator (`to_dnf` /
`to_dnf_neg` recognise `Element`/`NotElement` as dual literals over a container),
so `!Element[...]` and `NotElement[...]` normalise to the same form.

**Data structures.** Trivial: two `Expr` arguments copied into an `Element` node
that is handed to `evaluate`; the result is a single symbol or `NULL`. No
intermediate representation of its own.

**Complexity / limits.** As cheap as the underlying `Element` call. It decides
precisely when `Element` decides (e.g. `NotElement[I, Reals] -> True`,
`NotElement[2, Integers] -> False`) and is left unevaluated otherwise
(`NotElement[x, Reals]`); it performs no reasoning `Element` itself does not.
