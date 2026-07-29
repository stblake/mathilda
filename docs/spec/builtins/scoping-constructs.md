# Scoping Constructs

## Module
Implements lexical scoping by creating unique local variables.
- `Module[{x, y, ...}, expr]`
- `Module[{x = x0, ...}, expr]`

**Features**:
- `HoldAll`, `Protected`.
- Variables are renamed to `name$nnn` using `$ModuleNumber`.
- Created symbols have the `Temporary` attribute.

```mathematica
In[1]:= x = 1; Module[{x = 2}, x + 1]
Out[1]= 3

In[2]:= x
Out[2]= 1
```

## Block
Implements dynamic scoping by temporarily overriding symbol values.
- `Block[{x, y, ...}, expr]`
- `Block[{x = x0, ...}, expr]`

**Features**:
- `HoldAll`, `Protected`.
- Affects only values, not names.
- Restores original values and attributes after execution.

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

