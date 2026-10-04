---
source: src/piecewise.c
---
**Algorithm.** `builtin_ramp` is the positive part max(x, 0), the standard
spelling of a rectified linear unit. It classifies the sign with `ustep_class`:
a non-negative argument is returned unchanged, a negative one becomes zero, and
an undecidable or non-real argument is left unevaluated. The zero returned for a
negative argument carries the **argument's own exactness** — `Ramp[-1.]` is `0.`,
`Ramp[-3]` is the exact `0`, and an MPFR argument returns a zero at its precision —
so a Real vector maps to a Real vector and an integer vector to an integer one
with no mixed-head result (this is why `Ramp` needs no gate on its output where
`Clip` does).

**Data structures.** One `ustep_class` call and at most one `expr_copy`; no
intermediate expressions, unlike `UnitBox`. The ND kernel (`REG_U(Ramp)`) maps a
packed or visible numeric `NDArray` element-wise and preserves the element type,
so the buffer is answered in place.

**Complexity / limits.** `O(1)` per element; `Compile[]` lowers it at scalar and
rank-1 shapes. Like `UnitStep`, `Ramp` threads over an `Interval` argument (it is
non-decreasing, so endpoint threading is a rigorous enclosure). A genuinely
complex argument, or one whose sign cannot be certified, stays symbolic.
