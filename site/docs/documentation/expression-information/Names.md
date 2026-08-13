# Names

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Names["string"] gives a sorted list of the names of symbols matching the string. Names[patt] matches a string pattern with metacharacters * (zero or more characters) and @ (one or more non-uppercase characters), or a RegularExpression["re"]. Names[{p1, p2, ...}] matches any of the patterns. Names[] lists all symbol names.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Names["List*"]
Out[1]= {"List", "ListConvolve", "ListCorrelate", "ListPlot", "ListQ"}

In[2]:= Names["Ar@"]
Out[2]= {"Arg", "Array", "Arrow"}

In[3]:= Names[RegularExpression["Si."]]
Out[3]= {"Sin"}

In[4]:= MemberQ[Names["System`*"], "System`Sin"]
Out[4]= True
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

**Attributes:** `Protected`.

## References

**See also:** [RegularExpression](../../string-operations/RegularExpression/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
- Tests: [`tests/test_names.c`](https://github.com/stblake/mathilda/blob/main/tests/test_names.c)
