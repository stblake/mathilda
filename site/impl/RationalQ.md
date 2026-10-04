---
source: src/core.c
---
**Algorithm.** `builtin_rationalq` returns `True` when its one argument is an exact
rational number. That is exactly two cases, tested directly: the argument is
integer-like (`expr_is_integer_like` — an `EXPR_INTEGER` or a GMP `EXPR_BIGINT`,
since every integer is rational) or it has head `Rational` (`core_head_is(arg,
SYM_Rational)`, a `Rational[p, q]`). Otherwise `False`.

**Data structures.** None beyond the argument. Two cheap tests — a type tag check
and an interned-pointer head comparison — then a fresh `True`/`False` symbol.

**Complexity / limits.** `O(1)`. Reals and MPFR numbers are *not* rational here
(`2.5` is `False`) even when they have an exact rational value, because the test is
on representation; an irrational constant such as `Pi`, and any non-numeric
expression, is `False`. A call that is not single-argument declines and stays
symbolic.
