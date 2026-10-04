---
source: src/assoc_ops.c
---
**Algorithm.** `builtin_applyto` is the in-place update `x = f[x]`, returning the
new value. It is `HoldFirst`, so the first argument arrives as an unevaluated
l-value. `lvalue_root` walks it to the underlying symbol, accepting a bare symbol
`s`, a part `s[[i]]`, or an association entry `s[key]`; the builtin then
evaluates the l-value, applies `f`, and writes the result back by evaluating
`Set[lhs, newvalue]` — so the write-back works for all three l-value shapes
through the one `Set` path.

**Data structures.** Plain `Expr` trees; the current value is adopted into the
`f[...]` call, and the write-back is an ordinary `Set` evaluation. No hashing of
its own.

**Complexity / limits.** The cost is that of evaluating `f[x]` plus the
assignment. The l-value must already have a value — otherwise
`ApplyTo::rvalue` is emitted (via the message funnel) and nothing is changed;
a wrong argument count raises `ApplyTo::argrx`.
