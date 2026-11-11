# Scoping Constructs

Pattern substitution into a scoping construct is **capture-avoiding**. When a
rule's right-hand side holds a `Module`/`Block`/`With`/`Function`/`Table` and a
substituted value carries a free symbol whose name equals one of the construct's
bound locals, that local is alpha-renamed to a fresh symbol first, so the free
symbol is not swallowed: `g[v_] := Module[{e = 1}, v + e]; g[e + 1]` gives
`2 + e`. A bound local that shadows a pattern variable of the same name likewise
does not receive that variable: `sh[e_] := Module[{e = 1}, e + 1]; sh[99]` gives
`2`.

## Module
Implements lexical scoping by creating unique local variables.
- `Module[{x, y, ...}, expr]`
- `Module[{x = x0, ...}, expr]`

**Features**:
- `HoldAll`, `Protected`.
- Variables are renamed to `name$nnn` using `$ModuleNumber`.
- Created symbols have the `Temporary` attribute.
- **Temporary locals are removed once nothing refers to them** (v0.6056), as
  in the reference.
  - **Before.** Every local stayed for the session: about 550 bytes and one
    symbol per local per call, and +111,005 symbols (+44 MB) after one issue
    #112 Dirac trace.
  - **When.** A local is removed at its `Module`'s exit, or, if it escaped
    (returned, stored, captured by a closure, a rule or an UpValue), when the
    last outside reference is gone. That is checked as each top-level
    evaluation ends, and also during a long evaluation once a few thousand
    escaped locals are waiting (v0.6058). So a loop that stores and drops a
    `Module` result or closure per iteration runs in bounded memory: 200,000
    such iterations in one `Do` leave +0.9 MB, where they left +50 MB before.
  - **Nothing is left behind** (v0.6058). The removed symbol's name is freed
    too. Before, names were kept, at about 23 bytes per local.
  - **Exemptions.** A local whose `Temporary` attribute is cleared is kept, as
    is one created over a name that already existed.

```mathematica
In[1]:= x = 1; Module[{x = 2}, x + 1]
Out[1]= 3

In[2]:= x
Out[2]= 1
```

## Conditions at the end of a scoping construct
In a definition, a right-hand side that is `Module`, `With` or `Block`
whose body ends in `val /; cond` makes the rule **conditional**. "Ends in"
means through nested scopes and through the last statement of a
`CompoundExpression`.
- When the rule is tried, the body runs, and `cond` is evaluated inside the
  scope, so it sees the local variables.
- If `cond` is `True`, the rule's value is `val`. Otherwise the rule does not
  apply: the next rule is tried, or the call stays unevaluated. Side effects
  already made by the body remain.
- A later definition with the same left-hand side does not replace a
  conditional rule. Both are kept, in definition order.
- Outside a definition, `Module[{t = 1}, t /; t > 0]` is an ordinary
  `Condition` expression.

```mathematica
m[x_] := Module[{t = x^2}, t /; t > 3]
m[x_] := "fallback"
{m[1], m[2]}                                   (* {"fallback", 4} *)
```

## Block
Implements dynamic scoping by temporarily overriding symbol values.
- `Block[{x, y, ...}, expr]`
- `Block[{x = x0, ...}, expr]`

**Features**:
- `HoldAll`, `Protected`.
- Affects only values, not names.
- **Localises attributes and built-in rules too** (v0.6170), as in the reference: inside
  `Block[{f}, ...]` the symbol has no attributes and no C builtin, so `Attributes[f]` is `{}`,
  an `Orderless` `f[b, a]` stays unsorted, `Block[{Sin}, Sin[{0, 1}]]` is `Sin[{0, 1}]`, and
  `Block[{Message}, ...]` makes `Message` calls inert. Before, the attributes were saved and
  restored but never cleared. Inside `Block[{Times}, ...]` the reference still multiplies
  numbers (`Times[2, 3, b]` is `6 b`); Mathilda leaves a blocked `Times` entirely inert.
- Restores original values, attributes and builtins after execution.
- **The returned value is evaluated again once the symbols are restored**
  (v0.6053), as in the reference: `x = 5; Block[{x}, x + y]` is `5 + y`, and
  with `f[a_] := a^2`, `Block[{f}, f[3]]` is `9`. Before, the fixed-point stamp
  made while the symbol was blocked survived the restore, so these came back
  as `x + y` and `f[3]`.

## Internal`InheritedBlock
`Block` that keeps the current definitions. Implemented in `src/modular.c`
alongside `Block` (it shares the frame machinery, so `TimeConstrained` unwinds
and `Return` behave the same way).
- `Internal`InheritedBlock[{x, y = v, ...}, expr]`

**Features**:
- Each localised symbol starts with a *copy* of its current own values, down
  values, up values, attributes and builtin behaviour (`Block` starts from an
  unset symbol with no attributes). `y = v` then assigns a new local value.
- Everything the body does to the symbols (new definitions, `Unprotect`,
  `ClearAttributes`, ...) is undone on exit.
- The usual idiom: temporary rules for a system symbol,
  `Internal`InheritedBlock[{Plus}, Unprotect[Plus]; Plus[a, b] = 0; a + b]` is
  `0`, and afterwards `a + b` is unchanged and `Plus` is `Protected` again.
- Attributes: `HoldAll`, `Protected`.

```
f[x_] := x^2; y = 7;
Internal`InheritedBlock[{f, y}, f[x_] := x^3; y++; {f[2], y}]   (* {8, 8} *)
{f[2], y}                                                       (* {4, 7} *)
```

## With
Defines local constants by lexical substitution.
- `With[{x = x0, ...}, expr]`
- `With[{x := x0, ...}, expr]`

**Features**:
- `HoldAll`, `Protected`.
- Replaces occurrences of symbols in the body before evaluation.

```mathematica
In[1]:= x = 10; With[{x = 5}, x^2]
Out[1]= 25
```

## Unique
Generates fresh symbols that have never been used before.
- `Unique[x]` — a new symbol `x$nnn`, numbered by `$ModuleNumber` (the same
  counter and the same form `Module` uses for its locals).
- `Unique["x"]` — a new symbol `xnnn`, with no `$`.
- `Unique[]` — a new symbol `$nnn`.
- `Unique[{x, y, ...}]` — a list of fresh symbols, one per element.

**Features**:
- `Protected`.
- Every suffix is advanced until the generated name is unused. In a list,
  distinct sources share one number (`{a$5, b$5}`), while a repeated source gets
  a new number each time (`Unique[{a, a}]` is `{a$5, a$6}`); `$ModuleNumber` then
  moves past the largest number used.
- `Unique["x"]` and `Unique[]` draw from their own counter, not `$ModuleNumber`.
- A user assignment to `$ModuleNumber` is honoured, as for `Module`.
- An argument that is neither a symbol nor a string emits `Unique::usym` and is
  left unevaluated.
- Created symbols have the `Temporary` attribute.

```mathematica
In[1]:= {Unique[a], Unique[{a, a}], Unique[{a, b}], Unique["x"], Unique[]}
Out[1]= {a$1, {a$2, a$3}, {a$4, b$4}, x1, $2}
```

## Scoping and pure functions

`Module` / `Block` / `With` reach into the body of a nested pure function, so a
local may be read *and assigned* from inside a `&`:

```mathematica
In[1]:= Module[{c = 0}, Scan[(c = c + 1) &, Range[4]]; c]
Out[1]= 4

In[2]:= Module[{a = 2}, Map[(a #) &, {1, 2, 3}]]
Out[2]= {2, 4, 6}
```

A nested `Function` shadows an enclosing local of the same name, in both the
list and bare-symbol parameter forms — the parameter wins:

```mathematica
In[3]:= Module[{x = 5}, Function[x, x + 1][10]]
Out[3]= 11
```

Nested *named-parameter* `Function`s are not closures: an inner
`Function[e, ...]` does not capture an outer `Function`'s parameter. Inject it
with `With[]` when that is needed.

