---
source: src/match.c
---
**Algorithm.** `KeyValuePattern` is an inert pattern head consumed by the matcher,
not an evaluation builtin (it is interned and marked `Protected` in
`src/patterns.c`; the matching lives in `match_internal`). `KeyValuePattern[{k1 -> p1,
…}]` matches an association — or a list of rules — that contains, for each
requirement `ki -> pi`, an entry whose key matches `ki` and whose value matches `pi`.
The requirements are solved with backtracking (`kvp_match_reqs`), so a consistent
assignment is found whenever one exists even when requirements share bound variables
(e.g. `KeyValuePattern[{k_ -> _, _ -> k_}]`). An empty spec matches any association,
and a bare single rule is the one-requirement form.

**Data structures.** The requirement list is read directly off the pattern's
argument; bindings accumulate in the matcher's `MatchEnv`, with `env_rollback` undoing
a failed branch of the backtracking search.

**Complexity / limits.** Each requirement is matched against the subject's entries;
the backtracking is exponential only in the number of requirements that share
variables, which is small in practice. The subject must have head `Association` or
`List`, or the match fails immediately.
