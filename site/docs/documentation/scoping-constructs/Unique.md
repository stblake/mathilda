# Unique

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Unique[] generates a new symbol; Unique["x"] or Unique[x] uses a name prefix; Unique[{x, ...}] gives a list of fresh symbols. Each is Temporary and never previously used.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= {Unique[], Unique["x"], Unique[{a, b}]}
Out[1]= {$1, x2, {a3, b3}}
```

### Applications (4)

```mathematica
In[2]:= Unique[]
Out[2]= $1

In[3]:= Unique[x]
Out[3]= x2

In[4]:= Unique["c"]
Out[4]= c3

In[5]:= Unique[{p, q}]
Out[5]= {p4, q4}
```

## Implementation notes

**Algorithm.** `builtin_unique` (`src/modular.c`) generates fresh, never-before-used
symbols. `Unique[]` uses the prefix `"$"`; `Unique[x]` or `Unique["x"]` uses the
symbol's name or the string as prefix (`unique_prefix_of`); `Unique[{a, b, ...}]`
returns a list of fresh symbols, one per element, all sharing one numeric suffix
(the Wolfram behaviour). The suffix is drawn from the module counter
`module_number`, the same monotone source `Module` uses for its `x$n` temporaries.

`unique_make_batch` guarantees freshness by *scanning*: it advances
`module_number` until `prefix<n>` names no existing symbol for every prefix in
the batch, then takes that value and post-increments. Each generated symbol is
created with `expr_new_symbol` and tagged `ATTR_TEMPORARY`, and the user-visible
`$ModuleNumber` OwnValue is re-synced from the counter (`unique_sync_module_number`),
so `Unique` and `Module` share one numbering. An unusable prefix (not a symbol or
string) makes the batch fail and the head returns `NULL`.

**Attributes & limits.** `Unique` has no special attributes — its single argument
is evaluated normally. It accepts zero or one argument; the list form yields one
suffix-sharing symbol per element.

- `Protected`.
- The numeric suffix is drawn from the shared `$ModuleNumber` counter (the same
  source `Module` uses), advanced until the generated name is unused, so every
  result is genuinely fresh and distinct.
- Created symbols have the `Temporary` attribute.

**Attributes:** `Protected`.

## References

**See also:** [Module](../../scoping-constructs/Module/)

- Source: [`src/modular.c`](https://github.com/stblake/mathilda/blob/main/src/modular.c)
- Specification: [`docs/spec/builtins/scoping-constructs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/scoping-constructs.md)
- Tests: [`tests/test_series.c`](https://github.com/stblake/mathilda/blob/main/tests/test_series.c)

## Notes & additional examples

### Notes

`Unique[]` generates a brand-new symbol, `Unique[x]` or `Unique["x"]` uses the
given name as a prefix, and `Unique[{a, b, ...}]` returns a list of fresh symbols
sharing one numeric suffix. Each name is formed from a prefix and a number drawn
from the `$ModuleNumber` counter — the same counter `Module` uses for its
locals.

Freshness is guaranteed by construction: the counter is advanced until the
candidate name is unused for every prefix in the batch, so a returned symbol has
never named anything before. Each generated symbol is `Temporary`.
