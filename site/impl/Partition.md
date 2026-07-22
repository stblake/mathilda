---
source: src/list.c
---
**Algorithm.** `builtin_partition` splits a list into sublists of length `n` with offset `d`
(default `d = n`, i.e. non-overlapping blocks), via the recursive `partition_rec`. At each level
it reads the block size `n` and offset `d` for that level (a plain integer applies to level 0,
or a `List` gives a per-level spec), computes the number of full blocks `(len − n)/d + 1`, and
emits each sublist `args[i·d .. i·d + n)` wrapped in the list's head. An `UpTo[n]` size allows a
short final block. It recurses into each element so multi-level specs partition nested arrays.
Trailing partial blocks (when no `UpTo`) are dropped, following the no-padding default.
