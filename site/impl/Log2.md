---
source: src/logexp.c
---
**Algorithm.** `builtin_log2` is `fixed_base_log(res, "Log2", 2, log2)`, the exact
base-2 twin of `Log10`. Log2[z] = Log[2, z] = Log[z]/Log[2]. A positive finite
machine real goes straight to libm's `log2` — the one-ulp error of
`log(z)/log(2)` at exact powers of two (so `Log2[1024.]` would miss `10.`) is why
the dedicated libm call is used rather than a ratio of logarithms. Exact powers
of two, symbolic, negative, complex, arbitrary-precision, and zero arguments all
fall through to `Log[2, z]`, keeping the two spellings consistent.

**Data structures.** A `double` on the fast arm, else a two-argument `Log[2, z]`
tree. The ND kernel (`UK_FIXED_LOG(Log2, log2, ND_LN2)`, `REG_U`) runs a packed or
visible real `NDArray` element-wise through libm `log2`, escaping a negative or
complex element to `clog(z)/ln(2)`.

**Complexity / limits.** `O(1)` per element; `Compile[]` lowers it at scalar and
rank-1 shapes. `Log2[0]` is `-Infinity`; the negative axis returns the principal
complex value, both from `Log`.
