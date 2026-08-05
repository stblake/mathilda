# Control Flow

## Do
Evaluates an expression sequentially over an iteration range.
- `Do[expr, n]`: Evaluates `expr` `n` times.
- `Do[expr, {i, imax}]`: Evaluates `expr` with `i` from 1 to `imax`.
- `Do[expr, {i, imin, imax, di}]`: Evaluates `expr` with `i` taking values from `imin` to `imax` in steps of `di`.
- `Do[expr, {i, list}]`: Evaluates `expr` with `i` taking values from `list`.
- `Do[expr, spec1, spec2, ...]`: Evaluates `expr` looping over `spec1` internally.

**Features**:
- `HoldAll`, evaluating its body only after arguments are substituted.
- Employs exact dynamic iteration identical to `Table` but discards the evaluated results, returning `Null`.
- Supports explicit break states (`Return`, `Break`, `Continue`, `Throw`, `Abort`, `Quit`).
- Can execute an infinite loop using `Do[expr, Infinity]`.
- A body that is machine-numeric throughout takes an automatic fast path
  (`src/numloop.c`) that runs it as a double-stack program with no `Expr`
  allocation. The fast path is built without evaluating anything: loop-invariant
  terms are folded only when they are syntactically numeric (literals, `Pi`,
  `Sqrt[2]`, `p/q`, arithmetic and elementary heads). A body containing anything
  else — notably a call to a user-defined function — declines the fast path and
  runs interpreted, so `expr` is evaluated exactly as many times as the iterator
  specifies and its side effects fire exactly that often. The same rule governs
  `For`, `While`, `Nest`, `NestWhile`, `FixedPoint`, `Fold` and `Map`.

```mathematica
In[1]:= Do[Print[i], {i, 3}]
Out[1]= Null

In[2]:= Do[If[i == 3, Break[]]; Print[i], {i, 5}]
Out[2]= Null
```

## For
Executes a loop with an initialization, condition test, increment, and body.
- `For[start, test, incr, body]`: Evaluates `start`, then repeatedly evaluates `body` and `incr` until `test` fails to give `True`.
- `For[start, test, incr]`: Executes the loop with a `Null` body.

**Features**:
- Evaluates its arguments in a nonstandard way (sequence: `test`, `body`, `incr`).
- Has attribute `HoldAll`.
- `Break[]` exits the loop.
- `Continue[]` skips the rest of the body and proceeds to evaluating `incr`.
- Exits as soon as `test` fails.
- Returns `Null` unless an explicit `Return` is evaluated.

```mathematica
In[1]:= For[i=0, i<4, i++, Print[i]]
Out[1]= Null
```

## While
Evaluates a test expression and, while it yields `True`, repeatedly evaluates a body expression.
- `While[test, body]`: Evaluates `test`, then `body`, repeatedly, until `test` first fails to give `True`.
- `While[test]`: Executes the loop with a `Null` body. Useful when `test` itself has side-effects.

**Features**:
- Has attribute `HoldAll`; both `test` and `body` are re-evaluated each iteration.
- `Break[]` inside `body` exits the loop, yielding `Null`.
- `Continue[]` inside `body` skips the rest of `body` and returns to re-evaluating `test`.
- `Return[v]` inside `body` causes `While` to yield `v`.
- `Throw`, `Abort`, and `Quit` propagate unchanged.
- If the very first evaluation of `test` is not `True`, `body` is never evaluated.
- Returns `Null` unless an explicit `Return` is issued.

```mathematica
In[1]:= n = 1; While[n < 4, n = n + 1]; n
Out[1]= 4

In[2]:= {a, b} = {27, 6}; While[b != 0, {t1, t2} = {b, Mod[a, b]}; a = t1; b = t2]; a
Out[2]= 3

In[3]:= n = 1; While[True, If[n > 10, Break[]]; n = n + 1]; n
Out[3]= 11
```

## Break
Exits the nearest enclosing `Do`, `For`, or `While` loop.
- `Break[]`: Takes no arguments; the enclosing loop terminates and yields `Null`.

**Features**:
- Has attribute `Protected`.
- Takes effect as soon as it is evaluated (e.g. inside an `If` within the body),
  escaping only the *innermost* enclosing loop.
- Outside any loop, `Break[]` emits the message `Break::nofwd` and returns
  `Hold[Break[]]` (inert, so feeding it back does not re-trigger).

```mathematica
In[1]:= Do[Print[i]; If[i > 2, Break[]], {i, 10}]
1
2
3

In[2]:= For[i = 1, i <= 10, i++, If[i > 2, Break[]]]; i
Out[2]= 3
```

## Continue
Proceeds to the next iteration of the nearest enclosing `Do`, `For`, or `While` loop.
- `Continue[]`: Takes no arguments; the remainder of the current loop body is skipped.

**Features**:
- Has attribute `Protected`.
- Takes effect as soon as it is evaluated. In `Do` it advances the iterator and
  re-tests; in `For` it evaluates the increment step then re-tests; in `While` it
  re-evaluates the test.
- Outside any loop, `Continue[]` emits the message `Continue::nofwd` and returns
  `Hold[Continue[]]`.

```mathematica
In[1]:= r = 0; Do[If[EvenQ[i], Continue[]]; r += i, {i, 10}]; r
Out[1]= 25

In[2]:= r = 0; For[i = 1, i <= 10, i++, If[EvenQ[i], Continue[]]; r += i]; r
Out[2]= 25
```

## If
Evaluates condition and executes the corresponding branch.
- `If[condition, t, f]`: Gives `t` if `condition` evaluates to `True`, and `f` if it evaluates to `False`.
- `If[condition, t, f, u]`: Gives `u` if `condition` evaluates to neither `True` nor `False`.

**Features**:
- `HoldRest`, evaluating only the chosen branch.
- Remains unevaluated if the condition is undetermined and `u` is not provided.
- `If[condition, t]` returns `Null` if `condition` evaluates to `False`.

```mathematica
In[1]:= If[True, x, y]
Out[1]= x

In[2]:= If[a < b, 1, 0, Indeterminate]
Out[2]= Indeterminate
```

## Which
Selects a value based on the first satisfied test.
- `Which[test1, value1, test2, value2, ...]`: Evaluates each `test_i` in turn, returning the corresponding `value_i` for the first test that yields `True`.

**Features**:
- Has attribute `HoldAll`; tests and values are held until `Which` examines them.
- If every `test_i` evaluates to `False`, `Which` returns `Null`. `Which[]` (no arguments) likewise yields `Null`.
- If a `test_i` evaluates to something other than `True` or `False`, a `Which` containing that test (in evaluated form) plus the remaining elements is returned unevaluated.
- A trailing test of `True` acts as a default clause.
- An odd number of arguments is a usage error; the expression is returned unevaluated.

```mathematica
In[1]:= Which[False, a, True, b]
Out[1]= b

In[2]:= Which[1 < 0, a, x == 0, b, 0 < 1, c]
Out[2]= Which[x == 0, b, 0 < 1, c]

In[3]:= Which[a == 1, x, a == 2, b] /. a -> 2
Out[3]= b

In[4]:= sign[x_] := Which[x < 0, -1, x > 0, 1, True, Indeterminate]; sign /@ {-2, 0, 3}
Out[4]= {-1, Indeterminate, 1}
```

## Switch
Selects a value by matching an expression against a sequence of patterns.
- `Switch[expr, form_1, value_1, form_2, value_2, ...]`: Evaluates `expr`, then pattern-matches it against each `form_i` in turn, returning the `value_i` for the first match.

**Features**:
- Attribute `HoldRest`; the form/value pairs are held until `Switch` examines them.
- Each `form_i` is evaluated immediately before its match is tried; only the chosen `value_i` is evaluated.
- A trailing form of `_` (Blank) acts as a catch-all default clause.
- If no `form_i` matches `expr`, the call is returned unevaluated.
- Wrong arity (no form/value pair, or an odd number of arguments after `expr`) is a usage error; the expression is returned unevaluated.
- Pattern variables bound by `form_i` (e.g. `{x_, y_}`) are *not* substituted into `value_i`; the form acts purely as a discriminator.
- `Break`, `Return`, and `Throw` inside the chosen value propagate as they do in any other held context.

```mathematica
In[1]:= Switch[42, _Integer, "int", _Real, "real", _, "other"]
Out[1]= "int"

In[2]:= t[e_] := Switch[e, _Plus, Together, _Times, Apart, _, Identity]; t[(1+x)/(1-x) + x/(1+x)]
Out[2]= Together

In[3]:= Switch[#, 1, one, 2, two, _, other] & /@ {1, 2, 3}
Out[3]= {one, two, other}
```

## Piecewise
Represents a piecewise function defined by a list of `{value, condition}` clauses.
- `Piecewise[{{val_1, cond_1}, {val_2, cond_2}, ...}]`: Returns the `val_i` paired with the first `cond_i` that yields `True`.
- `Piecewise[{{val_1, cond_1}, ...}, default]`: Uses `default` if none of the `cond_i` apply.
- `Piecewise[conds]` automatically rewrites to `Piecewise[conds, 0]`.
- Attribute `HoldAll` — only the surviving `val_i` is evaluated (by the outer evaluator).

**Simplification semantics**:
- Conditions are evaluated left-to-right.
- `{val_i, False}` clauses are dropped.
- At the first `{val_i, True}` all later clauses (and the default) are dropped; the `True` clause becomes the unconditional final case.
- If all preceding conditions are literally `False`, the value at the first `True` is returned directly.
- Consecutive clauses with structurally equal values are merged: their conditions are combined with `Or`.

```mathematica
In[1]:= Piecewise[{{Sin[x]/x, x < 0}, {1, x == 0}}, -x^2/100 + 1]
Out[1]= Piecewise[{{Sin[x]/x, x < 0}, {1, x == 0}}, (-x^2)/100 + 1]

In[2]:= Piecewise[{{e1, True}, {e2, d2}, {e3, d3}}]
Out[2]= e1

In[3]:= Piecewise[{{a, d1}, {b, d2}, {c, False}, {d, d4}}, ef]
Out[3]= Piecewise[{{a, d1}, {b, d2}, {d, d4}}, ef]

In[4]:= Piecewise[{{a, d1}, {b, d2}, {b, d3}, {c, d4}}, ef]
Out[4]= Piecewise[{{a, d1}, {b, d2 || d3}, {c, d4}}, ef]

In[5]:= Piecewise[{{Sin[x]/x, x < 0}, {1, x == 0}}, -x^2/100 + 1] /. x -> 5
Out[5]= 3/4
```

## TrueQ
Tests whether an expression evaluates explicitly to `True`.
- `TrueQ[expr]`: Yields `True` if `expr` is `True`, and `False` otherwise.

## Boole
Indicator function (Iverson bracket) — converts logical values to 0/1.
- `Boole[expr]`: Yields `1` if `expr` is `True` and `0` if it is `False`. Stays unevaluated otherwise.

**Attributes**: `Listable`, `Protected`. Boole automatically threads over `List` arguments, so a vector of predicates becomes a vector of 0/1 indicators.

Useful for expressing integrals and sums over predicate regions, and for one-hot / dummy encoding of categorical variables in statistics. `Boole[expr]` is semantically equivalent to `If[expr, 1, 0]` but is `Listable` and so vectorises naturally.

```mathematica
In[1]:= {Boole[False], Boole[True]}
Out[1]= {0, 1}

In[2]:= Boole[{True, False, True, True, False}]
Out[2]= {1, 0, 1, 1, 0}

In[3]:= Boole[x]
Out[3]= Boole[x]

In[4]:= Total[Boole[# > 0 & /@ {-1, 2, -3, 4, 5}]]
Out[4]= 3
```

## ConditionalExpression
Symbolic construct representing a value that is defined only when a condition holds.
- `ConditionalExpression[expr, cond]`: Represents `expr` when `cond` is `True`.
- `ConditionalExpression[expr, True]`: Evaluates to `expr`.
- `ConditionalExpression[expr, False]`: Evaluates to `Undefined`.
- Nested forms collapse: `ConditionalExpression[ConditionalExpression[e, c1], c2]` reduces to `ConditionalExpression[e, c1 && c2]`.

**Attributes**: `Protected`.

```mathematica
In[1]:= ConditionalExpression[a, True]
Out[1]= a

In[2]:= ConditionalExpression[a, False]
Out[2]= Undefined

In[3]:= ConditionalExpression[x^2, x > 0]
Out[3]= ConditionalExpression[x^2, x > 0]

In[4]:= ConditionalExpression[ConditionalExpression[e, c1], c2]
Out[4]= ConditionalExpression[e, c1 && c2]
```

## Return
Exits the nearest enclosing scope or loop, yielding a value.
- `Return[expr]`: Yields `expr` from the innermost enclosing scope boundary.
- `Return[]`: Yields `Null` (shorthand for `Return[Null]`).
- `Return[expr, h]`: Yields `expr` from the *nearest enclosing boundary* whose head is the symbol `h`. Boundaries with a different head propagate the marker outward unchanged so that `h` can be reached.

**Recognised boundary heads**: `Function`, `Module`, `Block`, `With`, `Do`, `For`, `While`. `CompoundExpression` and other Hold-free heads (including `If`, `Which`, `Switch`) propagate `Return` through unchanged so it can reach the enclosing boundary.

**Features**:
- `Protected`. No Hold attributes — arguments are evaluated before the marker takes effect.
- `Return` takes effect as soon as it is evaluated, even when it appears inside other functions (Plus, Times, etc. in argument positions still see the substituted value, but at the top level the marker propagates immediately through `CompoundExpression`).
- 1-arg `Return[expr]` is consumed by the *innermost* boundary on the call stack. The 2-arg `Return[expr, h]` form lets the user skip past intervening boundaries to a specific enclosing construct.
- If no enclosing boundary matches, `Return[expr]` (or `Return[expr, h]`) survives at the top level as a literal expression.

```mathematica
In[1]:= Function[x, If[x > 0, Return[positive], Return[negative]]][3]
Out[1]= positive

In[2]:= Module[{}, Do[Return[5], {3}]]
Out[2]= 5

In[3]:= Module[{}, Do[Return[5, Module], {3}]]
Out[3]= 5

In[4]:= Module[{}, Do[Return[5, Block], {3}]]
Out[4]= Return[5, Block]

In[5]:= f[n_] := Module[{s = 0}, Do[s = s + i; If[s > 10, Return[i]], {i, 1, n}]]; f[10]
Out[5]= 5
```

## Catch / Throw
Non-local exit: `Throw` stops evaluation and hands a value to the nearest
enclosing `Catch`. Unlike `Return` (which only escapes a scope boundary), a
`Throw` propagates through *any* enclosing expression — `Plus`, `Times`,
function application, `Map`, `Sum`, `Table`, etc.

- `Throw[value]`: returns `value` to the nearest enclosing `Catch`.
- `Throw[value, tag]`: caught only by a `Catch[expr, form]` whose `form` matches `tag`.
- `Throw[value, tag, f]`: if uncaught, the top-level value is `f[value, tag]`.
- `Catch[expr]`: returns the argument of the first `Throw` generated while evaluating `expr`, or `expr` itself if none is thrown.
- `Catch[expr, form]`: returns `value` from the first `Throw[value, tag]` whose `tag` matches `form`; throws with a non-matching tag propagate to an outer `Catch`. `tag` is re-evaluated each time it is compared to `form`.
- `Catch[expr, form, f]`: returns `f[value, tag]` instead of `value`.

**Features**:
- `Throw` is `Protected`; `Catch` is `HoldFirst, Protected` (it drives evaluation of its body itself, so it can intercept a throw; `form` and `f` evaluate normally).
- Implemented by sentinel propagation through the evaluator's normal return paths (no `setjmp`/`longjmp`), so every frame runs its own cleanup — leak-free.
- The first `Throw` evaluated wins; a tagless `Throw[value]` is not caught by a form-`Catch`.
- An uncaught `Throw[value]`/`Throw[value, tag]` returns `Hold[Throw[...]]` with a `Throw::nocatch` message; an uncaught `Throw[value, tag, f]` returns `f[value, tag]`.

```mathematica
In[1]:= Catch[a; b; Throw[c]; d; e]
Out[1]= c

In[2]:= f[x_] := If[x > 10, Throw[overflow], x!]; Catch[f[2] + f[11]]
Out[2]= overflow

In[3]:= Catch[Do[If[i! > 10^10, Throw[i]], {i, 100}]]
Out[3]= 14

In[4]:= Catch[Throw[a, u], u]
Out[4]= a

In[5]:= Catch[Throw[v, tg], tg, {#1, #2} &]
Out[5]= {v, tg}
```

## Goto / Label
Imperative jumps within a `CompoundExpression`. `Goto[tag]` transfers control to
the `Label[tag]` in the compound expression the `Goto` appears in directly, then
in enclosing ones. `Label[tag]` marks a jump target; as a statement it evaluates
to `Null`.

- `Goto[tag]`: scans the enclosing `CompoundExpression`'s statements for
  `Label[tag]` and resumes evaluation there — a forward jump (skipping the
  statements in between) or a backward jump (forming a loop). If no matching
  `Label` is found in the current compound expression, the `Goto` propagates to
  the enclosing one; if none matches anywhere it is left unevaluated.
- `Label[tag]`: marks a point that `Goto[tag]` can jump to. It must appear as an
  explicit element of a `CompoundExpression`.

**Features**:
- Both are `Protected`. `tag` is evaluated (conventionally a literal symbol or
  integer) and compared structurally to each `Label`'s tag.
- Like `Catch`/`Throw`, `Goto` is implemented by sentinel propagation through the
  evaluator's normal return paths (no `setjmp`/`longjmp`), so a `Goto` fired
  inside a nested call (e.g. an `If` branch) still reaches the enclosing
  `CompoundExpression`. Leak-free.
- A `Goto` loop is a genuine loop with no artificial iteration cap; termination
  is the program's responsibility (as with `While`).
- A `Goto[tag]` that reaches the top level with no matching `Label` anywhere
  emits a `Goto::nolabel` message (stderr) and returns the inert `Goto[tag]`
  node. The message fires only when truly unmatched — a `Goto` that legitimately
  propagates from an inner to an outer `CompoundExpression` mid-evaluation is
  silent.

```mathematica
In[1]:= Module[{i = 0, s = 0}, Label[top]; i = i + 1; s = s + i;
          If[i < 5, Goto[top]]; s]
Out[1]= 15

In[2]:= f[a_] := Module[{x = 1., xp}, Label[begin];
          If[Abs[xp - x] < 10^-8, Goto[end]]; xp = x; x = (x + a/x)/2;
          Goto[begin]; Label[end]; x];
        f[2]
Out[2]= 1.41421
```


## Compile / CompiledFunction

`Compile[{arg, ...}, expr]` or `Compile[{{arg, _Type}, ...}, expr]` builds a
`CompiledFunction` object that evaluates `expr` over raw machine numbers,
bypassing the symbolic evaluator entirely. Argument types are `_Real` (the
default for a bare symbol), `_Integer`, or `_Complex`. Applying the object to
numeric arguments runs typed bytecode with no `Expr` allocation; a symbolic
argument, or a body outside the compilable subset, transparently falls back to
the interpreter, so a `CompiledFunction` always returns what the original
expression would.

- **Attributes:** `HoldAll`, `Protected`. The body is compiled in its raw
  (unevaluated) form; the argument symbols stay local, so a global value for one
  of them does not leak in.
- **Compilable subset** (shared with the internal engine behind NDSolve): full
  scalar arithmetic and comparisons, `Mod`/`Quotient`, integer/real/complex
  `Power`, all elementary functions and every special function that has a
  machine kernel (`Gamma`, `Erf`, `BesselJ`, `Zeta`, …), `If`, `Sum`/`Product`,
  `With`/`Module` locals with `Set`/`AddTo`/`TimesBy`/`Increment`/…,
  `Do`/`While`/`For`, `CompoundExpression`, the functional heads (see below),
  and machine arrays (see further below). Anything else (a user-defined
  function, exact symbolic algebra) routes that application through the
  interpreter fallback. Use `CompileDiagnostics` to find out which.
- **Functional heads.** `Nest`, `NestList`, `Fold`, `FoldList`, `FixedPoint`,
  `FixedPointList`, `NestWhile`, `NestWhileList`, `Map`, `Scan`, `Select`,
  `TakeWhile`, `LengthWhile`, `AllTrue`, `AnyTrue`, `NoneTrue`, `First`, `Last`
  and `Table` compile, as do the structural heads
  `Reverse`, `Sort`, `Accumulate`, `Flatten`, `Transpose`, `Take[a, n]` and
  `Drop[a, n]` — those last are delegated to the same `NDArray` entry points the
  interpreter uses, so their compiled subset is the interpreted one, and a spec
  outside it (a `Sort` comparator, `Take[v, {2, 4}]`) declines cleanly.

  The function argument may be `Function[u, body]`, `Function[{u, …}, body]`,
  `Function[body]` using `#`/`#1`/`#2`, a bare head (`Sin`, `Plus`, or any
  builtin with a machine kernel), a symbol holding another `CompiledFunction`,
  `Composition[…]`, or `Identity`. It must be one of those *at compile time* —
  there is no runtime function value — and its arity must match exactly, since
  a short call leaves a parameter symbolic and so is not a machine value. A
  `Slot` inside a **named** lambda (`Function[u, # + u]`) is not bound, matching
  the interpreter, which substitutes names only there.

  `Fold[f, x, v]`, `Fold[f, v]`, `Map[f, v]` and `Scan[f, v]` take a rank-1
  array. `FixedPoint[f, x]`, `FixedPoint[f, x, n]`, `FixedPoint[f, x, SameTest -> s]`
  and `NestWhile[f, x, test]` are all supported.

  These decline rather than answer differently:

  - `Table` needs **integer** iterators — a real iterator is walked by repeated
    addition in the interpreter, which a closed form does not reproduce
    exactly — and a **non-integer-valued body**, because a packed buffer has no
    integer dtype and `Table[i, {i, 1, n}]` holds exact `Integer`s.
    `Table[1.0 i, {i, 1, n}]` is the compilable spelling.
  - `Map` needs rank 1 (at rank ≥ 2 it maps over *rows*) and a result element
    type equal to the source's, because mapping over a packed array repacks with
    the source dtype.
  - `Fold` over an empty vector, and `Nest`/`NestList`/`FixedPoint` with a
    negative count, fall back: all are left unevaluated by the interpreter.
  - `NestList`, `FoldList`, `FixedPointList` and `NestWhileList` refuse an
    integer element type, for the same reason `Table` does — a packed buffer has
    no integer dtype.
  - `Select` and `TakeWhile` fall back when nothing is selected: an empty result
    has no packed form, so the interpreter answers with a `List` and a length-0
    array would not be the same value.
  - An unbounded `FixedPoint`/`NestWhile` that does not converge falls back at
    the same 10⁶-application cap the interpreter uses, and is then left
    unevaluated exactly as the interpreter leaves it.
  - `Do`, `While`, `For` and `Scan` answer `Null`, which the machine lattice has
    no room for, so they compile only where their value is discarded — inside a
    `CompoundExpression`, not as the whole body.
- **Machine arrays.** An argument spec `{v, _Real, 1}` (or `_Complex`, any rank)
  declares an array parameter. A `List` argument is packed into a flat machine
  buffer at the boundary and the result is unpacked back to a `List`; an
  `NDArray` argument is used in place and an `NDArray` comes back, so state kept
  packed across many calls converts once rather than per call. Elementwise
  expressions over whole arrays are fused into a single strip-mined pass.

  Inside a body:

  - **`Part` reads** accept the full spec vocabulary. One scalar subscript per
    axis (`v[[i]]`, `m[[i, j]]`, `t[[i, j, k]]`, negatives counting from the
    end, computed indices) lowers *inline* — no allocation, and each axis is
    range-checked separately, because `m[[1, ncols + 5]]` is inside the buffer
    but off the end of its row. Every other spec — `Span`, `All`, a list of
    positions, partial indexing like `m[[2]]` for a row, or any mixture such as
    `m[[k, 2 ;; 4]]` — is delegated to the same `ndarray_part` the interpreter
    uses, so the compiled answer *is* the interpreted one by construction, at
    the cost of allocating the result.
  - **`Part` assignment** with the same vocabulary: `u[[i, j]] = x`,
    `u[[2 ;; 4]] = 0.`, `u[[All, 2]] = 0.`, `u[[{1, 3}]] = 7.`, and
    `u[[1 ;; 3]] = w` from a matching array. `+=`, `-=`, `*=`, `/=` work on a
    scalar position. The target must be an array the program **owns** — a
    `Module`/`With` local — because the write goes into the buffer in place.
    Writing through an *argument* is not in the subset: argument arrays are
    borrowed, and for a `List` argument packed at the boundary the write would
    vanish silently. Copy it into a local first.
  - **`ConstantArray[v, n]`** and `ConstantArray[v, {n1, n2, …}]` create one,
    at any rank. The rank must be evident from the source; the dimensions are
    ordinary expressions evaluated per call. An *integer* fill is refused on
    purpose — `ConstantArray[0, n]` holds exact integer zeros in the interpreter
    and an `NDArray` has no integer dtype, so compiling it would answer
    differently, not just faster. Use `ConstantArray[0., n]`.
  - **`Module`/`With` locals may be arrays**, including as the result. A local
    initialised from an argument is a copy, matching the interpreter's value
    semantics.
  - Not in the subset: array-valued `If` branches, `Sum` accumulators or `Nest`
    state (each would duplicate a handle without duplicating ownership), and
    `Table` as an array constructor.
- **Association parameter bags (read-only).** An argument spec `{p, _Association}`
  (value element type defaults `_Real`) or `{p, _Association, _Real|_Integer|_Complex}`
  declares a read-only association parameter. The bag is borrowed — the program
  reads it and never mutates or frees it. Inside a body: `Lookup[p, key]`,
  `Lookup[p, key, default]`, `KeyExistsQ`/`KeyMemberQ`/`KeyFreeQ`, `Length[p]` and
  `Values[p]` (a packed vector). The default must be literal; a **constant-key**
  `Lookup` is O(1) and, being pure, is hoisted out of an enclosing loop and
  computed once. A **runtime-varying integer or real key** — the `Lookup[p, i]`-in-
  a-loop pattern — also compiles (O(1) per probe, no per-iteration allocation); a
  key past the machine-scalar boundary (a string, an array) makes the whole body
  fall back to the interpreter. A looked-up value that does not fit the
  declared type, or an absent key with no default, cleanly declines to the
  interpreter. A **constant** association (a literal `<|…|>`, or a global captured
  under auto-compilation) is folded at compile time, so a `Table`/`Plot` body that
  reads a captured global association compiles. Compiled code can also **produce**
  associations as first-class values: `KeyDrop[p, keys]`, `KeyTake[p, keys]` (keys
  literal) and `Counts[v]` (a machine array to an `element -> count` association)
  build a new association the compiled function can return, feed to any reader
  (`Lookup[KeyDrop[p, "x"], "y"]`, `Total[Values[Counts[v]]]`), or chain through
  another transform (`KeyTake[KeyDrop[p, a], b]`) to any depth. Higher-order
  transforms compile the embedded function into a callee run per value (no
  interpreter at runtime): `Map[f, assoc]` transforms each value (keys copied
  through), `Select[assoc, pred]` filters by value — `f`/`pred` may be a pure
  function (`#^2 &`), a `Function[…]`, or a bare head, and a `Map`/`Select` may
  itself consume a produced association (`Map[f, KeyDrop[p, k]]`). A callee
  outside the compilable subset makes the whole body fall back. `KeyUnion` (a key
  list), `PositionIndex` (list-valued entries) and the grouping family
  (`Merge`/`GroupBy`) are not in the subset and fall back to the interpreter.
- **Counted iterators** in `Do`/`Sum`/`Product` accept every integer-bounded
  spelling the interpreter does: `Do[body, n]` and `Do[body, {n}]` (repeat n
  times), `{i, hi}`, `{i, lo, hi}`, and `{i, lo, hi, di}` with a nonzero integer
  literal step. `Sum`/`Product` require a named iterator, matching the
  interpreter, which does not accept a bare count for them. A missing spelling
  is not merely a slower loop — it makes the whole surrounding body fall back to
  the interpreter.
- **Result** is a boxed machine number: an `Integer`, `Real`, or `Complex`
  (real part only if the imaginary part is zero), or `True`/`False` for a
  Boolean body.
- Applying with the wrong number of arguments leaves the application
  unevaluated. The object is reference-counted and leak-free.
- **Option `RuntimeAttributes`** — default `{}`; the only other setting is
  `Listable` (written either `RuntimeAttributes -> Listable` or
  `RuntimeAttributes -> {Listable}`). A `Listable` object threads over its
  list-valued arguments exactly as a `Listable` symbol does:

  - Threading happens over the levels a parameter does **not** consume, so for a
    scalar parameter any list threads, while a rank-*r* array parameter
    (`{v, _Real, r}`) takes an *r*-deep list whole and threads only over deeper
    ones — a `Listable Compile[{{v, _Real, 1}}, Total[v]]` maps over the rows of
    a matrix.
  - Nested lists thread level by level; `{}` threads to `{}`; a non-list
    argument is reused for every element; and lists of unequal length report
    `Thread::tdlen` and leave the application unevaluated.
  - The attribute belongs to the **object**, not to its bytecode: a body outside
    the compilable subset threads too, and each element independently takes
    either the compiled or the interpreted path (so a symbolic element in an
    otherwise numeric list still gives the symbolic answer for that element).
  - An `NDArray` argument threads by the same rank rule and the result is packed
    back into an `NDArray`, so a packed array and the `List` it packs give the
    same answer. A [packed `List`](packed-arrays.md) threads the same way and
    comes back packed, with the element type **re-sniffed** rather than forced to
    `Real` — so an integer-valued body over a packed integer list returns
    `Integer`s.
  - `Options[Compile]` reports the default and `SetOptions[Compile, …]` changes
    it. An unrecognised option, or a `RuntimeAttributes` setting other than the
    two above, leaves `Compile[…]` unevaluated rather than quietly ignoring it.
  - `CompilePrint` shows an `Attributes  Listable` line for such an object.

- **Machine integers.** Integer arithmetic is exact: every integer operation that
  can overflow an `int64` abandons the compiled call, and the interpreter re-runs
  the body and promotes to a bigint, so `Compile[{{n, _Integer}}, n^3][3000000]`
  gives `27000000000000000000` and not a wrapped value. Loop accumulators inherit
  this, as do the integer-CLOSED heads — `Power` with a non-negative exponent,
  `Factorial`, `Gamma`, `Binomial`, `Pochhammer`, `Fibonacci`, `LucasL`, `Im`,
  `Arg`, `FractionalPart` — each of which returns an **Integer**, matching the
  interpreter's head, and defers outside the domain where the interpreter returns
  an integer at all (`Factorial[-1]` is `ComplexInfinity`, `7^-3` is `1/343`).
  `GCD`, `LCM`, `PowerMod`, `Divisible`, `EvenQ`, `OddQ`, `IntegerLength` and
  `IntegerExponent` compile over integer arguments.

  `Sqrt`, `Divide` and the transcendentals still widen an integer argument to a
  Real, so `Compile[{{n, _Integer}}, n/2][3]` is `1.5` where the interpreter says
  `3/2`. That divergence is inherent rather than an oversight: Rationals and
  symbolic radicals are not machine numbers, and the Wolfram Language's `Compile`
  behaves the same way.

- **Packed array arguments.** A [packed `List`](packed-arrays.md) — an ordinary
  `List` that Mathilda stores as a dense buffer — is borrowed at the boundary and
  returned packed, so packing survives a compiled call:

  ```mathematica
  In[1]:= f = Compile[{{u, _Real, 1}}, u^2 + 1.];
          r = f[Range[1., 200000.]]; {NDArrayQ[r], Head[r]}
  Out[1]= {True, List}
  ```

  The result's presentation follows its argument: a plain `List` gives a plain
  `List`, a packed one gives a packed one (at any size — a *derived* array
  inherits presentation, so no threshold applies), and `NDArray[...]` gives
  `NDArray[...]`. A body that BUILDS its array (`ConstantArray`, `Table`,
  `NestList`) has nothing to inherit and follows the producer rule instead. A
  complex result is never packed, because `Complex[re, 0.]` is not a form the
  evaluator produces.

  A dtype mismatch — `Range[n]` infers `int64` and can reach a `_Real` parameter —
  costs one O(n) cast rather than declining the call to the interpreter. The cast
  declines rather than rounds: an `int64` magnitude past 2^53 has no exact
  `double`, and narrowing a float into an `_Integer` slot would change values.
  Both cases leave the interpreter to answer.

- **`$AutoCompilation`.** Many builtins compile a body *behind the caller's
  back*, purely as an optimisation: `Plot`, `Plot3D`, `ContourPlot`,
  `DensityPlot`, the parametric and vector plots, `Table` over an inexact
  iterator, `NIntegrate`, `NSum`, `FindRoot`, `NDSolve`, and the bodies of `Do`,
  `For`, `While`, `Map`, `Nest`, `Fold` and `FixedPoint`. `$AutoCompilation` is
  `True` by default and switches all of that off when set to `False`.

  ```mathematica
  In[1]:= $AutoCompilation = False; Table[i^2, {i, 1., 4.}]
  Out[1]= {1., 4., 9., 16.}
  ```

  An automatically compiled body is contracted to give the interpreter's answer,
  so this changes speed and nothing else. It is there so the two paths can be
  compared — a differential run flips it and diffs every output, and a user who
  suspects a compiled path has it wrong can confirm in one line.

  `Compile[]` and any `CompiledFunction` the user built are **not** affected:
  those were asked for. Only `True` or `False` is accepted; anything else is
  refused with an `$AutoCompilation::flagset` message and the symbol is rolled
  back to the live state, so reading it back never lies about which path is
  running. The environment variable `MATHILDA_NO_AUTOCOMPILE` starts a session
  with it off (and `MATHILDA_NO_NUMLOOP` covers the loop-body compiler alone).

  The companion switch for the other invisible optimisation is
  `$AutoArrayPacking` — see
  [`packed-arrays.md`](packed-arrays.md#turning-it-off--autoarraypacking).

- **`RuntimeOptions`.** `RuntimeOptions -> {"CatchMachineIntegerOverflow" ->
  False}` keeps the wrapped `int64` instead of deferring to the interpreter;
  `"Speed"` and `"Quality"` are shorthands for `False` and the default `True`.
  It is a **semantics** switch, not a performance one — the choice is compiled
  into each instruction rather than tested per call, so leaving checking on costs
  nothing measurable. Turning it off makes the object answer differently from the
  interpreter once a result leaves the machine-integer range, which is why it is
  opt-in. Auto-compiled paths (`Plot`, `Table`, `NIntegrate`, `NDSolve`) always
  check, and `SetOptions[Compile, …]` does not reach them.

- **Integer arrays.** `{v, _Integer, r}` is a rank-`r` machine-integer array
  parameter, alongside `_Real` and `_Complex`. `Range[n]`, an integer-valued
  `Table`, and `ConstantArray` with an integer fill build packed integer buffers
  and return Lists of **Integers**, matching the interpreter element for element
  and staying exact past 2^53. A `List` of Integers goes in and a `List` of
  Integers comes back — never an `NDArray`, since no user syntax builds an
  integer-typed one.

```mathematica
In[1]:= f = Compile[{{x, _Real}}, x^2 + 1]
Out[1]= CompiledFunction[{x}, x^2 + 1]

In[2]:= f[3.0]
Out[2]= 10.0

In[3]:= f[a]                 (* symbolic argument -> interpreter fallback *)
Out[3]= 1 + a^2

In[4]:= g = Compile[{{n, _Integer}}, Module[{s = 0.}, Do[s = s + 1/i^2, {i, 1, n}]; s]];
        g[100]
Out[4]= 1.63498

In[5]:= Compile[{{z, _Complex}}, z^2][1.0 + 2.0 I]
Out[5]= -3.0 + 4.0 I

In[6]:= Compile[{{m, _Real, 2}}, m[[All, 1]]][{{1., 2.}, {3., 4.}}]
Out[6]= {1.0, 3.0}

In[7]:= (* RuntimeAttributes -> Listable: the object threads over lists *)
        h = Compile[{{x, _Real}}, If[x > 0, 1., -1.], RuntimeAttributes -> Listable];
        h[{1., -2., 3.}]
Out[7]= {1.0, -1.0, 1.0}

In[8]:= (* a rank-1 parameter consumes one level, so this maps over the rows *)
        Compile[{{v, _Real, 1}}, Total[v], RuntimeAttributes -> Listable][
          {{1., 2.}, {3., 4.}}]
Out[8]= {3.0, 7.0}

In[9]:= (* a 5-point stencil: read an argument grid, write a local copy *)
        Compile[{{a, _Real, 2}},
          Module[{n = Length[a], b = a},
            Do[b[[i, j]] = (a[[i - 1, j]] + a[[i + 1, j]] +
                            a[[i, j - 1]] + a[[i, j + 1]])/4,
               {i, 2, n - 1}, {j, 2, n - 1}];
            b]][Table[1.0 (10 i + j), {i, 1, 3}, {j, 1, 3}]]
Out[9]= {{11.0, 12.0, 13.0}, {21.0, 22.0, 23.0}, {31.0, 32.0, 33.0}}
```

A worked end-to-end example — an explicit finite-difference solver for the 2-D
wave equation, verified against an exact discrete solution and benchmarked
against the interpreter, Wolfram Language and `NDSolve` — is summarised in
[`docs/design/compile_state.md`](../../design/compile_state.md).

## CompileDiagnostics

`CompileDiagnostics[argspec, expr]` reports whether `expr` compiles for the given
`Compile[]` argument specification, and if not, **which subexpression stopped
it**. `argspec` is exactly what `Compile` takes.

This exists because a bail is otherwise invisible. When the compiler meets a
construct it cannot lower it returns nothing, the caller quietly interprets, and
the answer is still correct — just an order of magnitude slower. And the cost is
not proportional: the compilable subset is a *cliff*, so a single unsupported
head anywhere in a body costs the **entire** body, including everything in it
that would have compiled.

- **Attributes:** `HoldAll`, `Protected`.
- **On success** the result carries `"Compiled" -> True`, the `"ResultType"`
  (`"Real"`, `"Integer"`, `"Complex"`, `"Boolean"`), the `"Instructions"` count,
  the number of `"CommonSubexpressions"` the optimiser hoisted, and
  `"InstructionsUnoptimized"` — the same body compiled with the optimiser off,
  so what code generation actually removed is visible.
- **On failure** it gives `"Compiled" -> False`, a `"Reason"`, and the
  `"Subexpression"` — the **innermost** node the emitter could not lower, which
  is the actual cause rather than the construct that happens to contain it.
- A malformed argspec leaves the call unevaluated, exactly as `Compile` does.

Setting the environment variable `MATHILDA_COMPILE_DIAG=1` makes the same report
appear on stderr whenever a numeric builtin (`Plot`, `NIntegrate`, `NSum`,
`ContourPlot`, …) falls back to the interpreter because its body did not compile.

```mathematica
In[1]:= CompileDiagnostics[{{x, _Real}}, Sin[x] + x^2]
Out[1]= {"Compiled" -> True, "ResultType" -> "Real", "Instructions" -> 4,
         "CommonSubexpressions" -> 0, "InstructionsUnoptimized" -> 4}

In[2]:= CompileDiagnostics[{{x, _Real}}, Sin[x] + BarnesG[x]]
Out[2]= {"Compiled" -> False,
         "Reason" -> "no machine lowering for this head at these argument types",
         "Subexpression" -> "BarnesG[x]"}

In[3]:= CompileDiagnostics[{{x, _Real}}, Sin[x] + y]
Out[3]= {"Compiled" -> False,
         "Reason" -> "symbol is not a declared argument and holds no machine value",
         "Subexpression" -> "y"}

In[4]:= CompileDiagnostics[{{z, _Complex}}, Zeta[z]]   (* real kernel, no complex one *)
Out[4]= {"Compiled" -> False,
         "Reason" -> "no machine lowering for this head at these argument types",
         "Subexpression" -> "Zeta[z]"}
```

## CompilePrint

`CompilePrint[cf]` prints the bytecode of a `CompiledFunction`. Where
`CompileDiagnostics` says how *much* code there is, this says **which** code —
the only way to see whether the optimiser folded the constants, whether an array
chain fused, or whether a map will actually fan out across cores.

- **Attributes:** `Protected`. Deliberately **not** `HoldAll`, unlike `Compile`
  and `CompileDiagnostics`: the argument has to evaluate down to the compiled
  object, so both `CompilePrint[Compile[…]]` and `f = Compile[…];
  CompilePrint[f]` work. Anything that is not a `CompiledFunction` is left
  unevaluated. Returns `Null`.
- **Header** — the signature, an `Attributes` line when the object carries a
  `RuntimeAttributes` setting (threading happens before any bytecode runs, so it
  must not be invisible), each argument's register and declared type (an argument
  the body never reads is marked `(unused)`), the result register and type, the
  sizes of the three register banks, and the instruction, CSE and parallel-loop
  counts.
- **Registers** are named by bank and numbered by frame slot: `R` scalar, `V`
  array handle, `T` strip-mining tile.
- **Each instruction** gives the opcode and its raw operands on the left and a
  readable rendering on the right. A `>` in the gutter marks a branch target.
- **No addresses appear.** Machine kernels are resolved back to their symbol
  names, and callee programs and parallel loops are numbered, so the output is
  stable enough to diff between two versions of a body.
- **An uncompiled object** has no bytecode to show, so it reports the bail
  reason and the offending subexpression instead — the same answer
  `CompileDiagnostics` gives about a body, asked about an object.

```mathematica
In[1]:= CompilePrint[Compile[{{x, _Real}}, x^2 + 2.5 x + 1]]
Signature   CompiledFunction[{x : Real}, x^2 + 2.5 x + 1]
Arguments   1
              R0   : Real         x
Result      R1 : Real
Registers   3 scalar, 0 array, 0 tile   (frame 3 slots)
Program     5 instructions, 0 CSE, all-Real fast path

    0  POWI_R     R1, R0, 2                       R1 = R0^2
    1  MUL_RK     R2, R0, 2.5                     R2 = R0 * 2.5
    2  ADD_R      R1, R1, R2                      R1 = R1 + R2
    3  ADD_RK     R1, R1, 1                       R1 = R1 + 1
    4  RET        R1                              return R1
```

`MUL_RK` and `ADD_RK` are the optimiser's `K_BINK` forms — the constant lives
*in* the instruction, so no register has to be materialised for it. Seeing
`CONST` instructions here instead would mean the folding pass did not engage.

An array body shows all three banks, the strip-mined loop and the fan-out
marker:

```mathematica
In[2]:= CompilePrint[Compile[{{v, _Real, 1}}, v^2 + 2 v + 1]]
Signature   CompiledFunction[{v : Real[1]}, v^2 + 2 v + 1]
Arguments   1
              V0   : Real[1]      v
Result      V4 : Real[1]
Registers   4 scalar, 1 array, 3 tile   (frame 200 slots)
Program     19 instructions, 0 CSE, 1 parallel loop

    0  A_SIZE     R1, V0                          R1 = Length[V0] (flat)
    1  A_NEWLIKE  V4, V0, V0  [a:arr b:arr -> Real]  V4 = buffer like V0 (Real)
    2  CONST      R2, 0                           R2 = 0
    3  GT_IK      R3, R1, 0                       R3 = R1 > 0
    4  JZ         R3, -> 18                       if !R3 goto 18
    5  APAR       R2, R1, -> 18  <ploop #0>       parallel R2 over [0, R1), else fall through; join 18
>   6  VSETLEN    R2, R1                          vlen = min(R1 - R2, 64)
    7  VLOAD_R    T5, V0, R2                      T5 = V0[R2 ...]
    ...
   17  LOOP       R2, R1, -> 6                    if ++R2 < R1 goto 6
>  18  RET        V4                              return V4
```

Special functions lower to machine kernels, which are named rather than printed
as pointers:

```mathematica
In[3]:= CompilePrint[Compile[{{x, _Real}}, Gamma[x^2 + 2.5]]]
    ...
    2  KERN_RR    R1, R1, <Gamma>                 R1 = Gamma[R1]
    3  RET        R1                              return R1

In[4]:= CompilePrint[Compile[{x}, Integrate[x, x]]]
Signature   CompiledFunction[{x : Real}, Integrate[x, x]]
Program     not compiled — every call runs the interpreter
Reason      no machine lowering for this head at these argument types
Bailed on   Integrate[x, x]
```
