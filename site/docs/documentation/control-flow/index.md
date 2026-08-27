# Control Flow

29 built-in function(s) in this category.

- [`$Epilog`]($Epilog.md) — $Epilog  _(Stable)_
- [`$Post`]($Post.md) — $Post  _(Stable)_
- [`$Pre`]($Pre.md) — $Pre  _(Stable)_
- [`$PrePrint`]($PrePrint.md) — $PrePrint  _(Stable)_
- [`$PreRead`]($PreRead.md) — $PreRead  _(Stable)_
- [`$RecursionLimit`]($RecursionLimit.md) — $RecursionLimit  _(Stable)_
- [`Boole`](Boole.md) — Boole[expr]  _(Stable)_
- [`Break`](Break.md) — Break[] exits the nearest enclosing Do, For, or While loop.  _(Stable)_
- [`Catch`](Catch.md) — Catch[expr]  _(Stable)_
- [`Compile`](Compile.md) — Compile[{x, ...}, expr] or Compile[{{x, _Real}, ...}, expr] builds a CompiledFunction that evaluates expr over machine numbers (types _Real, _Integer, _Complex; default _Real), falling back to the interpreter for symbolic arguments or non-compilable bodies. With RuntimeAttributes -> Listable the object threads over List arguments; the default is RuntimeAttributes -> {}. RuntimeOptions -> {"CatchMachineIntegerOverflow" -> False} (or the shorthand RuntimeOptions -> "Speed") lets machine-integer arithmetic wrap instead of falling back to the interpreter, which is faster and gives a different answer from the interpreter once a result leaves the machine-integer range; the default True never does. WorkingPrecision -> n compiles real/complex arithmetic in MPFR at n decimal digits (one fixed precision for the whole function), for the straight-line arithmetic + elementary-function subset; MachinePrecision (the default) keeps the machine path unchanged. "BigIntegers" -> True makes integer arithmetic exact (GMP) instead of int64.  _(Stable)_
- [`CompileDiagnostics`](CompileDiagnostics.md) — CompileDiagnostics[argspec, expr] reports whether expr compiles for the given Compile[] argument specification, and if not, the innermost subexpression that could not be lowered. Accepts the same WorkingPrecision -> n / "BigIntegers" -> True options as Compile[], so it also reports whether the arbitrary-precision subset lowers (ResultType MPFRReal/MPFRComplex/BigInteger). For a compiled body it also gives the result type and the instruction count with and without the optimiser.  _(Stable)_
- [`CompilePrint`](CompilePrint.md) — CompilePrint[cf] prints the bytecode of the CompiledFunction cf: its argument and result registers with their types, the scalar/array/tile register banks, and one line per instruction giving both the raw operands and a readable rendering. For an object whose body did not compile it reports the bail reason instead. Returns Null.  _(Stable)_
- [`ConditionalExpression`](ConditionalExpression.md) — ConditionalExpression[expr, cond]  _(Stable)_
- [`Continue`](Continue.md) — Continue[] proceeds to the next iteration of the nearest enclosing Do, For, or While loop.  _(Stable)_
- [`Do`](Do.md) — Do[expr, n] evaluates expr n times.  _(Stable)_
- [`Equivalent`](Equivalent.md) — Equivalent[e1, e2, ...]  _(Stable)_
- [`For`](For.md) — For[start, test, incr, body] executes start, then repeatedly evaluates body and incr until test fails to give True.  _(Stable)_
- [`Goto`](Goto.md) — Goto[tag]  _(Stable)_
- [`If`](If.md) — If[cond, t]  _(Stable)_
- [`Implies`](Implies.md) — Implies[p, q]  _(Stable)_
- [`Label`](Label.md) — Label[tag]  _(Stable)_
- [`Piecewise`](Piecewise.md) — Piecewise[{{val_1, cond_1}, {val_2, cond_2}, ...}]  _(Stable)_
- [`Return`](Return.md) — Return[expr]  _(Stable)_
- [`Switch`](Switch.md) — Switch[expr, form_1, value_1, form_2, value_2, ...]  _(Stable)_
- [`Throw`](Throw.md) — Throw[value]  _(Stable)_
- [`TrueQ`](TrueQ.md) — TrueQ[expr] yields True if expr is True, and False otherwise.  _(Stable)_
- [`Which`](Which.md) — Which[test1, value1, test2, value2, ...]  _(Stable)_
- [`While`](While.md) — While[test, body] evaluates test, then body, repeatedly, until test first fails to give True.  _(Stable)_
- [`Xor`](Xor.md) — Xor[e1, e2, ...]  _(Stable)_
