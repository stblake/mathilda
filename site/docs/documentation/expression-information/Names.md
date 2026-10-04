# Names

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Names["string"] gives a sorted list of the names of symbols matching the string. Names[patt] matches a string pattern with metacharacters * (zero or more characters) and @ (one or more non-uppercase characters), or a RegularExpression["re"]. Names[{p1, p2, ...}] matches any of the patterns. Names[] lists all symbol names.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Names["List*"]
Out[1]= {"List", "ListConvolve", "ListCorrelate", "ListGradient", "ListInterpolation", "ListPlot", "ListQ"}

In[2]:= Names["Ar@"]
Out[2]= {"Area", "Arg", "Array", "Arrow"}

In[3]:= Names[RegularExpression["Si."]]
Out[3]= {"Sin"}

In[4]:= MemberQ[Names["System`*"], "System`Sin"]
Out[4]= True
```

### Applications (4)

```mathematica
In[5]:= Names["Sin"]
Out[5]= {"Sin"}
```

@ matches a run of one-or-more non-uppercase letters

```mathematica
In[6]:= Names["Ar@"]
Out[6]= {"Area", "Arg", "Array", "Arrow"}
```

A PCRE pattern, anchored to the whole name

```mathematica
In[7]:= Names[RegularExpression["Si."]]
Out[7]= {"Sin"}
```

A backtick pattern enumerates builtins

```mathematica
In[8]:= MemberQ[Names["System`*"], "System`Sin"]
Out[8]= True
```

## Algorithm

names.c - Names[] and friends: enumerate symbol-table names by pattern.

```text
  Names["string"]            names matching a string pattern
  Names[patt]                names matching an arbitrary string pattern patt
  Names[{p1, p2, ...}]       names matching any of the p_i
  Names[]                    all names in the symbol table
```

A string pattern is matched against the whole name (anchored) and supports two metacharacters:

```text
  *   matches zero or more characters
  @   matches one or more characters that are NOT uppercase letters
```

Every other character (including the ` used in context prefixes) is literal.

A pattern element may instead be RegularExpression["re"], matched against the

```text
whole name via the PCRE2 engine (src/strings/regex).  When the engine is
```

unavailable the call emits Names::regavail and stays unevaluated.

The result is a List of Strings sorted ascending by byte value (strcmp),

```text
matching Wolfram-Language ordering for the common ASCII case.  All symbols in
```

the table are candidates -- there is no filtering of internal helper symbols.

Context handling: symbols are stored under bare names for the System` and Global` contexts (builtins and unqualified user symbols) and under explicit

```text
backtick-qualified names for other contexts.  A pattern element that itself
```

contains a backtick (e.g. "System`*") is matched against -- and returns -- each symbol's fully context-qualified name (System`Sin, Global`x, ...); a plain pattern (no backtick) is matched against, and returns, the stored name

```text
exactly as before.  This is what makes Names["System`*"] enumerate the
```

builtins instead of returning {} (nothing is stored with a literal "System`" prefix).

## Implementation notes

**Algorithm.** `builtin_names` (`src/names.c`) enumerates the symbol table by
pattern. Each argument element is compiled by `build_pat` into either a glob
(`EXPR_STRING`) or a PCRE2 program (`RegularExpression["re"]`, anchored as
`\A(?:re)\z`); a `List` argument becomes a set of alternative patterns, and no
argument means "match every name". The glob matcher `wl_glob_match` is a recursive
backtracking whole-string match with two metacharacters — `*` (zero or more of
anything) and `@` (one or more non-uppercase characters) — every other byte,
including the context `` ` ``, literal. `symtab_for_each` then visits every symbol,
`match_one` decides the emitted string, and the collected names are `qsort`ed by
`expr_compare` so that `Names[p]` is identical to `Sort[Names[p]]`.

**Data structures.** A `NameCollect` growable `char**` buffer accumulates owned
name strings during the visitor callback; each becomes an `EXPR_STRING` leaf of the
returned `List`. Context qualification is lazy: `full_qualified_name` prefixes a
bare name with `System`` when its `SymbolDef` has a `builtin_func` (or is a
kernel-interned System symbol) and `Global`` otherwise, caching the result across
the pattern loop for one symbol.

**Complexity / limits.** `O(S · P · L)` for `S` symbols, `P` patterns, and glob
cost `L` per name (regex cost dominates when a `RegularExpression` is used). A
backtick in a glob switches both the match target and the returned string to the
fully-qualified name, which is what lets `Names["System`*"]` list the builtins. A
`RegularExpression` pattern on a build without PCRE2 emits `Names::regavail`
(through `mth_message`) and leaves the call unevaluated. Attributes `Protected`.

**Attributes:** `Protected`.

## References

**See also:** [RegularExpression](../../string-operations/RegularExpression/)

- Source: [`src/names.c`](https://github.com/stblake/mathilda/blob/main/src/names.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
- Tests: [`tests/test_names.c`](https://github.com/stblake/mathilda/blob/main/tests/test_names.c)

## Notes & additional examples

### Notes

`Names["patt"]` gives the canonically sorted list of symbol names matching a
string pattern. The pattern is matched against the **whole** name (anchored) with
two metacharacters: `*` for zero or more of anything, and `@` for one or more
characters that are not uppercase letters (so `"Ar@"` stops before the capital in
`ArcSin`). Every other character, the context `` ` `` included, is literal; a
pattern element may instead be `RegularExpression["re"]`.

A pattern containing a backtick is matched against, and returns, each symbol's
fully context-qualified name (`System`Sin`, `Global`x`), which is what makes
`Names["System`*"]` enumerate all builtins. The result is always sorted, so
`Names[patt]` equals `Sort[Names[patt]]`. On a build without PCRE2 a
`RegularExpression` pattern emits `Names::regavail` and stays unevaluated.
