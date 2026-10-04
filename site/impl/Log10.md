---
source: src/logexp.c
---
**Algorithm.** `builtin_log10` is `fixed_base_log(res, "Log10", 10, log10)`.
Log10[z] is defined, as in Mathematica, as Log[10, z] = Log[z]/Log[10]. One fast
arm handles a positive finite machine real directly with libm's `log10`, because
the quotient of two separately rounded logarithms is wrong by an ulp at exact
powers — `log(1000.)/log(10.)` is 2.9999999999999996 where `log10(1000.)` is
exactly 3. Every other argument — an exact power of ten (detected exactly), a
symbolic z, a negative or complex argument, an arbitrary-precision number, or
zero — is forwarded to `Log[10, z]`, so the two spellings can never disagree.

**Data structures.** A single `double` on the fast arm; otherwise a two-argument
`Log[10, z]` expression handed back to the evaluator. The ND kernel
(`UK_FIXED_LOG(Log10, log10, ND_LN10)`, registered with `REG_U`) maps a packed or
visible real `NDArray` element-wise through libm `log10`; a negative or complex
element escapes to `clog(z)/ln(10)`, so exact powers of ten stay exact on the
buffer as well.

**Complexity / limits.** `O(1)` per element; `Compile[]` lowers it at scalar and
rank-1 shapes. `Log10[0]` is `-Infinity` and the negative real axis returns the
principal complex value, both inherited from `Log`.
