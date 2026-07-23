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

