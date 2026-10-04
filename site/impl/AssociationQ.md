---
source: src/assoc.c
---
**Algorithm.** `builtin_associationq` takes a single argument and returns the
symbol `True` or `False` according to `assoc_is_wellformed` (`assoc_struct.h`).
Having the `Association` head is not sufficient: the constructor leaves a
malformed `Association[1, 2]`, or a `Map`/`Apply` result such as
`Association[f[a -> 1]]`, unevaluated, and such a node is an ordinary expression,
not an association. The well-formedness test therefore requires every entry to be
a two-argument `Rule` or `RuleDelayed`.

**Data structures.** A pure structural scan over the argument's child array; no
hash index is built and the argument is only read, never copied.

**Complexity / limits.** `O(n)` in the number of entries (a bare `<|...|>`
literal that already carries a validated shape still re-scans). A bare list of
rules `{a -> 1}` answers `False` — it is a `List`, not an `Association`. The
head is a `*Q` predicate: it always returns a Boolean and never stays
unevaluated.
