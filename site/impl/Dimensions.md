---
source: src/core.c
---
**Algorithm.** `builtin_dimensions` (in `src/core.c`) measures the shape of a rectangular nested structure all of whose levels share the same head (taken from the top-level expression's head). The recursive helper `get_dimensions` records each level's `arg_count`, then recurses into the first child to get the candidate sub-shape and verifies every sibling has identical depth and dimensions; as soon as the structure becomes ragged it stops and returns the dimensions found so far. An optional second argument caps the depth (`Infinity` maps to the internal cap).

**Data structures.** Fixed-size `int64_t dims[DIMENSIONS_MAX_DEPTH]` stack buffers (`DIMENSIONS_MAX_DEPTH = 64`) hold per-level extents; the result is a `List` of integers.
