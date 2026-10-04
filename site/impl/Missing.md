---
source: src/assoc_ops.c
---
**Algorithm.** `Missing` is an inert, `Protected` head that *represents* missing data;
it has no evaluation builtin, so `Missing[]`, `Missing["reason"]` and
`Missing["reason", data]` persist unchanged. It is produced by the association
machinery: `Lookup` and key access return `Missing["KeyAbsent", key]` for an absent
key (`make_missing` builds it), `KeyUnion` fills padded gaps with it, and
`JoinAcross` fills unmatched rows with `Missing["Unmatched"]`.

**Data structures.** A plain `Expr` function node with head `Missing`; the second
argument, when present, carries the key or reason that explains the absence.

**Complexity / limits.** None — it is a symbolic marker. Test for it with `MissingQ`
and drop such entries from an association with `DeleteMissing`; `Lookup`'s
third-argument default is the usual way to supply a value in place of one.
