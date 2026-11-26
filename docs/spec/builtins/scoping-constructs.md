# Scoping Constructs

Pattern substitution into a scoping construct is **capture-avoiding**. When a
rule's right-hand side holds a `Module`/`With`/`Function` and a substituted
value carries a free symbol whose name equals one of the construct's bound
locals, that local is alpha-renamed to a fresh symbol first, so the free symbol
is not swallowed: `g[v_] := Module[{e = 1}, v + e]; g[e + 1]` gives `2 + e`. A
bound local that shadows a pattern variable of the same name likewise does not
receive that variable: `sh[e_] := Module[{e = 1}, e + 1]; sh[99]` gives `2`.

`Block` and the iterator heads `Table`, `Do`, `Sum` and `Product` are **not**
lexical: they localise dynamically, so a substituted value is placed into their
body and iterator specs unchanged, and the iterator then binds it, as in the
reference: `g[f_] := Sum[f, {k, 1, 3}]; g[1/k]` gives `11/6`, and
`h[f_] := Table[f, {k, 3}]; h[k^2]` gives `{1, 4, 9}`. For the same reason an
enclosing `With` or `Function` substitutes into an iterator spec too:
`With[{k = 2}, Hold[Table[k, {k, 3}]]]` is `Hold[Table[2, {2, 3}]]`.

## Module
Implements lexical scoping by creating unique local variables.
- `Module[{x, y, ...}, expr]`
- `Module[{x = x0, ...}, expr]`

**Features**:
- `HoldAll`, `Protected`.
- Variables are renamed to `name$nnn` using `$ModuleNumber`.
- **`$ModuleNumber` is a counter, read through `Internal`ModuleNumber[]`** (v0.6653).
  The symbol holds the delayed OwnValue `$ModuleNumber :> Internal`ModuleNumber[]`, so
  `Module` and `Unique` advance the counter without rewriting a definition or moving the
  evaluation clock. A user `$ModuleNumber = n`, or `Clear[$ModuleNumber]`, takes effect
  at the next `Module` or `Unique`. A fresh local's initial value is bound without a
  clock bump.
  - Divergences: `OwnValues[$ModuleNumber]` shows the sentinel rule. Modules run inside
    `Block[{$ModuleNumber = n}, …]` leave the counter where they left it.
  - See `docs/design/scoping_substitution.md`.
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
- **A malformed local is rejected** (audit Block#2): an entry that is neither a symbol nor `sym = v` (`Module[{f[1]}, 2]`) issues `Module::lvsym` and the call stays unevaluated. Local names of any length stay unique (audit Module#1).

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
- **A duplicate local is rejected** (audit Block#1): `Block[{x = 1, x = 2}, x]` issues
  `Block::dup` and stays unevaluated, leaving the outer `x` untouched; `Module::dup` and
  `With::dup` likewise. Before, the frame saved `x` twice and its exit restored `x = 1`.
- **The returned value is evaluated again once the symbols are restored**
  (v0.6053), as in the reference: `x = 5; Block[{x}, x + y]` is `5 + y`, and
  with `f[a_] := a^2`, `Block[{f}, f[3]]` is `9`. Before, the fixed-point stamp
  made while the symbol was blocked survived the restore, so these came back
  as `x + y` and `f[3]`.
- **A malformed local is rejected** (audit Block#2): `Block[{1}, 2]` issues `Block::lvsym` and stays unevaluated.

## With

Capture-avoiding: a nested `Function`, `Module` or `With` binding a name that occurs free in a substituted value has that local renamed first: `With[{y = x}, Function[x, x + y][10]]` is `10 + x` (audit With#1).
Defines local constants by lexical substitution.
- `With[{x = x0, ...}, expr]`
- `With[{x := x0, ...}, expr]`

**Features**:
- `HoldAll`, `Protected`.
- Replaces occurrences of symbols in the body before evaluation.
- **Rules and definitions scope their pattern variables** (audit With#2): a pattern
  variable named like a local is renamed to `x$` throughout that rule instead of having the
  value substituted into the pattern name, so `With[{x = 5}, x_ -> x]` is `x$_ -> x$` and
  `With[{x = 5}, f[x_] := x]; f[2]` is `2` (Module renames the same way). Before, these gave
  the invalid pattern `Pattern[5, _] -> 5` and `5`.

```mathematica
In[1]:= x = 10; With[{x = 5}, x^2]
Out[1]= 25
```
- **A malformed local is rejected** (audit With#3): a bare symbol, `1 = 2` or `f[1] = 2` issues `With::lvset` and the call stays unevaluated.
- **Several lists nest** (audit With#4): `With[{x = 1}, {y = x + 1}, y]` is `2`.

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
- `Unique["str"]` needs a valid symbol name: `Unique["1abc"]` issues `Unique::usym` (audit Unique#1).
- `Unique[x, attrs]` (and `Unique[{x, ...}, attrs]`) gives the fresh symbols the attributes `attrs` as well as `Temporary` (audit Unique#2).

- The `$ModuleNumber` counter is not bounded by the machine range: past
  `2^63 - 1` it keeps counting in GMP, so `$ModuleNumber = 2^63 - 1;
  {Unique[w], Unique[v]}` is `{w$9223372036854775807, w$9223372036854775808}`
  (the second was `v$-9223372036854775808`), and `$ModuleNumber = 10^30` is
  honoured by `Unique` and `Module` alike.

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

