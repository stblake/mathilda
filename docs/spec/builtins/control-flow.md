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

## TrueQ
Tests whether an expression evaluates explicitly to `True`.
- `TrueQ[expr]`: Yields `True` if `expr` is `True`, and `False` otherwise.

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

