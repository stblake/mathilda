# Expression Information

## Attributes
Returns the core evaluation attributes assigned to a symbol.
- `Attributes[symbol]`

**Features**: 
- Common attributes include `Flat` (associativity), `Orderless` (commutativity), `Listable` (automatic threading over lists), `HoldFirst`, `HoldRest`, `HoldAll` (evaluation control), and `Protected`.

```mathematica
In[1]:= Attributes[Plus]
Out[1]= {Flat, Listable, NumericFunction, OneIdentity, Orderless}
```

## ClearAttributes
Removes attributes from symbols.
- `ClearAttributes[s, attr]`: Removes `attr` from the list of attributes of the symbol `s`.
- `ClearAttributes["s", attr]`: Removes `attr` from the attributes of the symbol named `"s"` if it exists.
- `ClearAttributes[s, {attr1, attr2, ...}]`: Removes several attributes at a time.
- `ClearAttributes[{s1, s2, ...}, attrs]`: Removes attributes from several symbols at a time.

**Features**:
- `HoldFirst`, `Protected`.
- `ClearAttributes` modifies `Attributes[s]`.
- Cannot clear attributes of a `Locked` symbol.
- Clearing an attribute that is not set is a no-op.

```mathematica
In[1]:= SetAttributes[f, Listable]
In[2]:= f[{1, 2, 3}]
Out[2]= {f[1], f[2], f[3]}
In[3]:= ClearAttributes[f, Listable]
In[4]:= f[{1, 2, 3}]
Out[4]= f[{1, 2, 3}]

In[5]:= SetAttributes[f, {Flat, Orderless, OneIdentity}]
In[6]:= ClearAttributes[f, OneIdentity]
In[7]:= Attributes[f]
Out[7]= {Flat, Orderless}

In[8]:= ClearAttributes[f, {Flat, Orderless}]
In[9]:= Attributes[f]
Out[9]= {}

In[10]:= SetAttributes[{g, h}, Protected]
In[11]:= ClearAttributes[{g, h}, Protected]
In[12]:= Attributes[g]
Out[12]= {}
```

## AtomQ, NumberQ, IntegerQ, MachineNumberQ
Predicates for testing expression types.
- `AtomQ[expr]`: `True` if the expression has no parts.
- `NumberQ[expr]`: `True` if the expression is a numeric type (Integer, Real, Rational, Complex).
- `IntegerQ[expr]`: `True` if the expression is an Integer.
- `MachineNumberQ[expr]`: `True` if `expr` is a machine-precision (IEEE
  double) real, or a `Complex` of two finite machine-precision reals.
  Returns `False` for exact numbers (Integer, BigInt, Rational), for
  arbitrary-precision reals (`MPFR` / `1.0`50`), and for non-finite
  doubles (`Infinity`, `Indeterminate`, or an `EXPR_REAL` carrying
  `inf`/`NaN`). Attributes: `Protected`.

```mathematica
In[1]:= MachineNumberQ[Sin[1000.]]
Out[1]= True

In[2]:= MachineNumberQ[Exp[1000.]]      (* overflows to +inf *)
Out[2]= False

In[3]:= MachineNumberQ[-29037945.290347]
Out[3]= True

In[4]:= MachineNumberQ[N[Pi, 30]]       (* MPFR, not machine *)
Out[4]= False

In[5]:= MachineNumberQ[1.0 + 2.0 I]
Out[5]= True

In[6]:= MachineNumberQ[1 + 2 I]         (* exact Gaussian integer *)
Out[6]= False
```

## PossibleZeroQ
Hybrid symbolic-numeric test for whether `expr` is identically zero.
The general problem is undecidable (Richardson 1968); `PossibleZeroQ`
is fast but not always exact — it may produce a small number of
one-sided false positives on identities that hold only on a measure-
zero subset of the complex plane.

Pipeline (each stage exits early on a definite verdict):

1. **Structural** — literal `0`, `Complex[0, 0]`, named non-zero
   constants like `Pi` are decided in O(1).
2. **Rational normalisation** — `Together`, `Cancel`, `Expand` plus
   `is_zero_poly` decides every identity in `Q(x_1, …, x_n)` exactly.
3. **Numeric precision ladder** — for closed-form numeric inputs,
   numericalize at machine precision and compare against a
   cancellation-aware threshold. Bump to MPFR (200/500/1000 bits) and
   re-evaluate when the residual lies near the rounding-noise floor.
4. **Schwartz–Zippel** — for inputs with free symbols, substitute
   random real/complex samples and recurse into the numeric stage.
   Requires four independent confirmations.

Attributes: `Listable`, `Protected`.

```mathematica
In[1]:= PossibleZeroQ[E^(I Pi/4) - (-1)^(1/4)]
Out[1]= True

In[2]:= PossibleZeroQ[(x + 1)(x - 1) - x^2 + 1]
Out[2]= True

In[3]:= PossibleZeroQ[(E + Pi)^2 - E^2 - Pi^2 - 2 E Pi]
Out[3]= True

In[4]:= PossibleZeroQ[E^Pi - Pi^E]
Out[4]= False

In[5]:= PossibleZeroQ[2^(2 I) - 2^(-2 I) - 2 I Sin[Log[4]]]
Out[5]= True

In[6]:= PossibleZeroQ[Sqrt[x^2] - x]
Out[6]= False    (* fails for negative x or complex x off the cut *)

In[7]:= PossibleZeroQ[Sin[x]^2 + Cos[x]^2 - 1]
Out[7]= True
```

## $MachinePrecision, $MachineEpsilon, $MinMachineNumber, $MaxMachineNumber, $MaxNumber, $MinNumber
Read-only system constants describing the floating-point range and
granularity Mathilda was compiled against. Each is a `Protected` symbol
holding an `OwnValue`; attempts to redefine them emit
`Set::wrsym`.

- `$MachinePrecision`: number of decimal digits of precision in a
  machine-precision number. Derived from the platform's `DBL_MANT_DIG`
  (typically `53 Log[10, 2] ~ 15.9546` on IEEE 754).
- `$MachineEpsilon`: the gap between `1.0` and the next-larger
  machine-precision number. Equals `DBL_EPSILON`.
- `$MinMachineNumber`: the smallest normalized positive machine-precision
  number. Equals `DBL_MIN`.
- `$MaxMachineNumber`: the largest finite machine-precision number.
  Equals `DBL_MAX`.
- `$MaxNumber` / `$MinNumber`: the largest / smallest positive value
  representable at machine precision in the arbitrary-precision backend.
  In `USE_MPFR` builds these reflect MPFR's current emax/emin one step
  in from infinity / zero; without MPFR they collapse onto
  `$MaxMachineNumber` / `$MinMachineNumber`.

```mathematica
In[1]:= $MachinePrecision
Out[1]= 15.9546

In[2]:= $MachineEpsilon
Out[2]= 2.22045*10^-16

In[3]:= {$MinMachineNumber, $MaxMachineNumber}
Out[3]= {2.22507*10^-308, 1.79769*10^308}

In[4]:= MachineNumberQ[$MaxNumber]   (* MPFR, not machine *)
Out[4]= False
```

## ListQ, VectorQ, MatrixQ
Predicates for testing lists and their structures.
- `ListQ[expr]`: `True` if the head of `expr` is `List`.
- `VectorQ[expr]`: `True` if `expr` is a list, none of whose elements are themselves lists.
- `VectorQ[expr, test]`: `True` if `expr` is a list and `test` yields `True` for each element.
- `MatrixQ[expr]`: `True` if `expr` is a list of lists of the same length, with no deeper lists.
- `MatrixQ[expr, test]`: `True` if `expr` is a matrix and `test` yields `True` for each element.

## EvenQ, OddQ
- `EvenQ[n]`: `True` if `n` is an even integer.
- `OddQ[n]`: `True` if `n` is an odd integer.

## Identity
- `Identity[expr]` returns `expr` unchanged. Useful as a default callable in higher-order functions.
- Attributes: `Protected`.

```mathematica
In[1]:= Identity[x]
Out[1]= x

In[2]:= Identity[1 + 2]
Out[2]= 3

In[3]:= Map[Identity, {a, b, c}]
Out[3]= {a, b, c}
```

## Composition
- `Composition[f1, f2, f3, ...]` represents the symbolic composition of the functions `f1, f2, f3, ...`. Applied to arguments, the composition acts innermost-first:
  `Composition[f, g, h][x, y]` -> `f[g[h[x, y]]]`.
- Attributes: `Flat`, `OneIdentity`, `Protected`.
- Algebraic simplifications applied automatically:
  - `Composition[]` -> `Identity`
  - `Composition[f]` -> `f`
  - `Identity` arguments are dropped: `Composition[f, Identity, g]` -> `Composition[f, g]`
  - Adjacent inverse pairs cancel: `Composition[f, InverseFunction[f]]` -> `Identity`, and similarly for the reverse order.
- Infix syntax: `f1 @* f2 @* ...` parses as `Composition[f1, f2, ...]` (precedence 625, left-associative). The `Flat` attribute then flattens any nested `Composition` calls produced by the parser.

```mathematica
In[1]:= Composition[f, g, h][x, y]
Out[1]= f[g[h[x, y]]]

In[2]:= f @* g @* h @ x
Out[2]= f[g[h[x]]]

In[3]:= Composition[f, Identity, g]
Out[3]= Composition[f, g]

In[4]:= Composition[f, InverseFunction[f]][x]
Out[4]= x

In[5]:= Composition[f, g] @* Composition[a, b]
Out[5]= Composition[f, g, a, b]
```

## ComposeList
- `ComposeList[{f1, f2, ..., fn}, x]` generates the list `{x, f1[x], f2[f1[x]], ..., fn[...f2[f1[x]]...]}`. The output has `n + 1` elements.
- Attributes: `Protected`.
- ComposeList builds the symbolic applications and lets the evaluator reduce them in the normal way, so it composes naturally with pure functions and explicit operator-style sequences.

```mathematica
In[1]:= ComposeList[{a, b, c, d}, x]
Out[1]= {x, a[x], b[a[x]], c[b[a[x]]], d[c[b[a[x]]]]}

In[2]:= ComposeList[{f, g}[[{1, 2, 1, 1, 2}]], x]
Out[2]= {x, f[x], g[f[x]], f[g[f[x]]], f[f[g[f[x]]]], g[f[f[g[f[x]]]]]}

In[3]:= ComposeList[{1 - # &, 1/# &}[[{2, 2, 1, 2, 2, 1}]], x]
Out[3]= {x, 1/x, x, 1 - x, 1/(1 - x), 1 - x, x}
```

## Print
- `Print[expr1, expr2, ...]` prints the expressions to stdout and returns `Null`.
- Supports `FullForm` and `InputForm` wrappers.

```mathematica
In[1]:= Print["Result: ", x + y]
Result: x + y
Out[1]= Null

In[2]:= Print[x + y // FullForm]
Plus[x, y]
Out[2]= Null
```

## FullForm
- `FullForm[expr]` causes `expr` to be printed in its full internal form.
- `expr // FullForm` is a common shortcut.


## HoldForm
- `HoldForm[expr]` prints as the expression `expr`, with `expr` maintained in an unevaluated form.
- Has attribute `HoldAll` and `Protected`.
- It dissolves prior to printing but preserves the unevaluated `expr` it contains.

```mathematica
In[1]:= HoldForm[1 + 1]
Out[1]= 1 + 1
```
## Evaluate
Causes `expr` to be evaluated even if it appears as the argument of a function whose attributes specify that it should be held unevaluated.
- `Evaluate[expr]`

**Features**:
- `Protected`.
- `Evaluate` only overrides `HoldFirst`, `HoldRest`, and `HoldAll` attributes when it appears directly as the head of the function argument that would otherwise be held.
- `Evaluate` does not override `HoldAllComplete`.
- `Evaluate` works only on the first level, directly inside a held function. It does not penetrate into deeper subexpressions.
- Outside of held contexts, `Evaluate` acts as identity.

```mathematica
In[1]:= Evaluate[1+1]
Out[1]= 2

In[2]:= Hold[Evaluate[1+1], 2+2]
Out[2]= Hold[2, 2 + 2]

In[3]:= Hold[Evaluate[1+1], Evaluate[2+2], Evaluate[3+3]]
Out[3]= Hold[2, 4, 6]

In[4]:= Hold[f[Evaluate[1+2]]]
Out[4]= Hold[f[Evaluate[1 + 2]]]

In[5]:= Hold[Evaluate[Sin[Pi/6]]]
Out[5]= Hold[1/2]

In[6]:= Hold[Evaluate[2^10]]
Out[6]= Hold[1024]

In[7]:= HoldForm[Evaluate[1+1]]
Out[7]= 2

In[8]:= Hold[Evaluate[Length[{a,b,c}]]]
Out[8]= Hold[3]

In[9]:= Evaluate[Head[{1,2,3}]]
Out[9]= List
```

## ReleaseHold
Removes `Hold`, `HoldForm`, `HoldPattern`, and `HoldComplete` in `expr`.
- `ReleaseHold[expr]`

**Features**:
- `Protected`.
- `ReleaseHold` removes only one layer of `Hold` etc.; it does not remove inner occurrences in nested `Hold` etc. functions.
- `ReleaseHold` traverses into subexpressions of `expr` and strips any hold wrapper it finds, but does not recurse into the contents of the stripped wrapper.
- When `expr` does not contain any hold wrappers, `ReleaseHold` acts as identity.
- `ReleaseHold` removes all standard unevaluated containers: `Hold`, `HoldForm`, `HoldComplete`, and `HoldPattern`.

```mathematica
In[1]:= Hold[1+1]
Out[1]= Hold[1 + 1]

In[2]:= ReleaseHold[Hold[1+1]]
Out[2]= 2

In[3]:= ReleaseHold /@ {Hold[1+2], HoldForm[2+3], HoldComplete[3+4]}
Out[3]= {3, 5, 7}

In[4]:= ReleaseHold[f[Hold[1+2]]]
Out[4]= f[3]

In[5]:= ReleaseHold[f[Hold[1+g[Hold[2+3]]]]]
Out[5]= f[1 + g[Hold[2 + 3]]]

In[6]:= ReleaseHold[Hold[Hold[1+1]]]
Out[6]= Hold[1 + 1]

In[7]:= ReleaseHold[Hold[Sin[Pi/6]]]
Out[7]= 1/2

In[8]:= ReleaseHold[{f[Hold[1+2]], g[HoldForm[3+4]]}]
Out[8]= {f[3], g[7]}
```

## HoldPattern
Keeps a pattern expression in unevaluated form while still allowing it to act as a pattern for matching.
- `HoldPattern[expr]`

**Features**:
- Attributes: `{HoldAll, Protected}`.
- `HoldPattern[p]` is equivalent to `p` in the pattern matcher; the matcher transparently unwraps a single-argument `HoldPattern` before matching.
- Useful on the left-hand side of rules and assignments, because those positions are normally evaluated before being used for matching. Wrapping in `HoldPattern` stops that evaluation and preserves the literal pattern shape.
- `HoldPattern` is removed by one layer of `ReleaseHold`.

```mathematica
In[1]:= HoldPattern[_+_] -> 0
Out[1]= HoldPattern[Blank[] + Blank[]] -> 0

In[2]:= a + b /. HoldPattern[_+_] -> 0
Out[2]= 0

In[3]:= Cases[{a -> b, c -> d}, HoldPattern[a -> _]]
Out[3]= {a -> b}

In[4]:= MatchQ[a + b, HoldPattern[_+_]]
Out[4]= True
```

## Unevaluated
Represents the unevaluated form of `expr` when it appears as the argument to a function.
- `Unevaluated[expr]`

**Features**:
- Attributes: `{HoldAllComplete, Protected}`.
- `f[Unevaluated[expr]]` passes `expr` to `f` as if `f` temporarily held that single argument; the `Unevaluated` wrapper is then stripped before `f`'s body runs, effectively yielding `f[expr]` with `expr` unevaluated.
- The wrapper is **not** stripped when the enclosing function holds the argument (e.g. `f` has `HoldAll`, or `HoldFirst`/`HoldRest` applies to that position).
- The wrapper is **not** stripped when the enclosing function has `HoldAllComplete`.
- Stripping happens **after** `Sequence` flattening, so a `Sequence` directly inside `Unevaluated` survives into the argument slot (`Length[Unevaluated[Sequence[a, b]]]` gives `2`).
- Nested `Unevaluated` wrappers are stripped one layer per evaluation step.
- As a top-level expression, `Unevaluated[expr]` evaluates to itself (because of `HoldAllComplete`).

```mathematica
In[1]:= Length[Unevaluated[Plus[5, 6, 7, 8]]]
Out[1]= 4

In[2]:= Length[Unevaluated[Sequence[a, b]]]
Out[2]= 2

In[3]:= Hold[Evaluate[Unevaluated[1+2]]]
Out[3]= Hold[Unevaluated[1 + 2]]

In[4]:= SetAttributes[f, HoldAll]; f[Unevaluated[1+2]]
Out[4]= f[Unevaluated[1 + 2]]

In[5]:= HoldComplete[Unevaluated[1+2]]
Out[5]= HoldComplete[Unevaluated[1 + 2]]

In[6]:= Attributes[Unevaluated]
Out[6]= {HoldAllComplete, Protected}
```

## HoldComplete
Shields `expr` completely from the standard evaluation process.
- `HoldComplete[expr1, expr2, ...]`

**Features**:
- Attributes: `{HoldAllComplete, Protected}`.
- `HoldComplete` prevents argument evaluation, `Sequence` flattening inside its own arguments, `Unevaluated` wrapper stripping, and `Evaluate` from firing. `Evaluate` cannot override `HoldAllComplete`.
- Structural substitution (via `ReplaceAll`, `Replace`, `ReplacePart`, etc.) still descends into `HoldComplete` because substitution is not part of evaluation.
- `HoldComplete` is removed by one layer of `ReleaseHold`.
- `HoldComplete` is a milder form of `Unevaluated` at top level: `HoldComplete` always keeps the wrapper, while `Unevaluated` is typically stripped by the enclosing function.

```mathematica
In[1]:= Attributes[HoldComplete]
Out[1]= {HoldAllComplete, Protected}

In[2]:= HoldComplete[1+1, Evaluate[1+2], Sequence[3, 4]]
Out[2]= HoldComplete[1 + 1, Evaluate[1 + 2], Sequence[3, 4]]

In[3]:= HoldComplete[Sequence[a, b]]
Out[3]= HoldComplete[Sequence[a, b]]

In[4]:= HoldComplete[f[1+2]] /. f[x_] :> g[x]
Out[4]= HoldComplete[g[1 + 2]]

In[5]:= ReleaseHold[HoldComplete[Sequence[1, 2]]]
Out[5]= Sequence[1, 2]
```

## HoldAllComplete (attribute)
An attribute that specifies that **all** arguments to a function are not to be modified or looked at in any way during evaluation. It is stricter than `HoldAll`.

A function with `HoldAllComplete`:
- does not evaluate any argument,
- does not flatten `Sequence[...]` that appears inside an argument,
- does not strip `Unevaluated[...]` wrappers inside its arguments,
- is not overridden by `Evaluate[...]` wrappers inside its arguments,
- does not apply any up-values associated with its arguments.

`HoldComplete` and `Unevaluated` are the two standard built-in heads that carry `HoldAllComplete`.

## InputForm
- `InputForm[expr]` causes `expr` to be printed in a form suitable for input (standard form in Mathilda).

## ToString
- `ToString[expr]`: returns a `String` containing the printed form of `expr` in `InputForm`.
- `ToString[expr, form]`: returns the printed form for the specified form. Supported forms are `InputForm` (default), `FullForm`, and `TeXForm`. `StandardForm` and `OutputForm` are accepted as aliases for `InputForm`.

**Features**:
- `Protected`.
- An unsupported form leaves the call unevaluated (e.g. `ToString[x, FooForm]` returns `ToString[x, FooForm]`), so a typo is visible at the call site rather than silently downgraded.

```mathematica
In[1]:= ToString[x^2 + y^3]
Out[1]= "x^2 + y^3"

In[2]:= ToString[x^2 + y^3, FullForm]
Out[2]= "Plus[Power[x, 2], Power[y, 3]]"

In[3]:= ToString[x^2 + y^3, TeXForm]
Out[3]= "x^{2}+y^{3}"
```

## ToExpression
- `ToExpression[input]`: parses the string `input` as Mathilda input (`InputForm`) and returns the resulting expression after evaluation.
- `ToExpression[input, form]`: uses interpretation rules corresponding to the specified form. Currently `InputForm`, `FullForm`, and `StandardForm` are accepted (all route through the same parser, which is form-agnostic).
- `ToExpression[input, form, h]`: wraps the head `h` around the parsed expression before it is returned to the evaluator. Using `h = Hold` produces an unevaluated `Hold[...]` wrapper.

**Features**:
- `Protected`, `Listable`. `ToExpression[{"1+1", "2+2"}]` evaluates to `{2, 4}`.
- Returns the symbol `$Failed` if the parser cannot consume the input.
- A non-string input or an unsupported `form` leaves the call unevaluated.

```mathematica
In[1]:= ToExpression["1+1"]
Out[1]= 2

In[2]:= ToExpression["1+1", InputForm, Hold]
Out[2]= Hold[1 + 1]

In[3]:= ToExpression["x+"]
Out[3]= $Failed
```

## Symbol
Refers to a symbol with the specified name, creating it if necessary.
- `Symbol["name"]`: returns the symbol named `"name"`.

**Features**:
- `Protected`.
- Every expression's `Head` matches `Symbol` for symbols; `x_Symbol` patterns therefore match any symbol.
- The string must satisfy the standard symbol-name syntax: each segment (separated by backticks) starts with a letter or `$`, followed by letters, digits, or `$`.
- A leading backtick (`Symbol["\`x"]`) makes the name relative to the current `$Context`. An embedded backtick (`Symbol["a\`x"]`) is treated as an absolutely-qualified name. A bare name is resolved through the standard `$Context` / `$ContextPath` rules.
- Invalid names emit `Symbol::symname` to `stderr` and leave the call unevaluated; non-string arguments also leave the call unevaluated.

```mathematica
In[1]:= Symbol["x"]
Out[1]= x

In[2]:= Head[%]
Out[2]= Symbol

In[3]:= {f[x], f["x"], f[2]} /. f[s_Symbol] :> g[s]
Out[3]= {g[x], f["x"], f[2]}

In[4]:= Symbol["a`x"]
Out[4]= a`x

In[5]:= Symbol["1x"]
Symbol::symname: The string "1x" cannot be used for a symbol name. A symbol name must start with a letter followed by letters and numbers.
Out[5]= Symbol["1x"]
```

## Information
Returns a formatted string containing the syntax and description of a symbol.
- `Information[symbol]`
- `?symbol` (shortcut)

```mathematica
In[1]:= ?Range
Out[1]= "Range[n]
	generates the list {1, 2, 3, ..., n}.
Range[n, m]
	generates the list {n, n + 1, ..., m - 1, m}.
Range[n, m, d]
	uses step d."
```

## MemberQ
- `MemberQ[list, form]`: Returns `True` if an element of list matches form, and `False` otherwise.
- `MemberQ[list, form, levelspec]`: Tests all parts of list specified by levelspec.
- `MemberQ[form]`: Represents an operator form of `MemberQ` that can be applied to an expression.

**Features**:
- `Protected`.
- Default option: `Heads -> False`.
- `form` can be a structural pattern.
- The first argument of `MemberQ` can have any head, not necessarily `List`.
- Returns immediately upon finding the first match.
- Standard level specifications are supported. The default value for `levelspec` in `MemberQ` is `{1}`.

```mathematica
In[1]:= MemberQ[{1, 3, 4, 1, 2}, 2]
Out[1]= True

In[2]:= MemberQ[{x^2, y^2, x^3}, x^_]
Out[2]= True

In[3]:= MemberQ[{{1, 1, 3, 0}, {2, 1, 2, 2}}, 0, 2]
Out[3]= True

In[4]:= MemberQ[{{1, 1, 3, 0}, {2, 1, 2, 2}}, 0]
Out[4]= False
```

## FreeQ
Yields `True` if no subexpression in `expr` matches `form`, and yields `False` otherwise.
- `FreeQ[expr, form]`
- `FreeQ[expr, form, levelspec]`

**Features**:
- `Protected`.
- By default, explores levels `{0, Infinity}` and option `Heads -> True` is enabled.
- `form` can be a structural pattern.

```mathematica
In[1]:= FreeQ[{1, 2, 4, 1, 0}, 0]
Out[1]= False

In[2]:= FreeQ[{a, b, b, a, a, a}, _Integer]
Out[2]= True

In[3]:= f[c_ x_, x_] := c f[x, x] /; FreeQ[c, x]
In[4]:= {f[3 x, x], f[a x, x], f[(1 + x) x, x]}
Out[4]= {3 f[x, x], a f[x, x], f[x (1 + x), x]}
```

## LeafCount
Gives the total number of indivisible subexpressions in an expression.
- `LeafCount[expr]`

**Features**:
- `Protected`.
- Counts the number of subexpressions in `expr` that correspond to "leaves" on the expression tree.
- By default `Heads -> True` includes the head of expressions and their parts. With `Heads -> False`, it excludes them.
- Evaluates atoms like `Rational` and `Complex` based on their structural representation as functions.

```mathematica
In[1]:= LeafCount[1 + a + b^2]
Out[1]= 6

In[2]:= LeafCount[f[a, b][x, y]]
Out[2]= 5

In[3]:= LeafCount[{1/2, 1 + I}]
Out[3]= 7
```

## ByteCount
Gives the number of bytes used internally by Mathilda to store the expression.
- `ByteCount[expr]`

**Features**:
- `Protected`.
- Uses `sizeof()` in C and measures the internal AST memory allocation boundaries, dynamically capturing sizes of individual strings, symbols, allocated blocks, arrays, and expression structs.

