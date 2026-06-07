---
source: src/patterns.c
---
`builtin_memberq` (`src/patterns.c`) returns `True`/`False` for whether any part of the first argument matches the second (a pattern), tested by `do_member_at_level` which applies the pattern matcher `match` at each position within the level spec (default level 1 — immediate elements). It supports integer/`{min,max}`/`All`/`Infinity` level specs and the `Heads` option; the one-argument form returns an operator `Function[MemberQ[#1, patt]]`.
