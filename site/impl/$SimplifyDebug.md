---
source: src/simp/simp_util.c
---
A diagnostic flag (not a builtin), given an OwnValue defaulting to `False` in
`simp_init`. When set to `True`, `simp_debug_enabled` (read directly off the
OwnValue list to avoid re-evaluation) causes `traced_call_unary` /
`simp_debug_log` to emit one stderr line per transform invocation inside the
Simplify search, in the form `/<TransformName>/: <input> -> <output> [<ms> ms]`,
used to diagnose slow Simplify calls and runaway candidate explosion.
