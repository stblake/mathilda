# StringEndsQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringEndsQ["string", patt]`**

Gives True if a suffix of "string" matches the string expression patt, and False otherwise.

**`StringEndsQ["string", {p1, p2, ...}]`**

Gives True if a suffix matches any of the pi.

**`StringEndsQ[{s1, s2, ...}, patt]`**

Gives the list of results for each of the si.

**`StringEndsQ[patt]`**

Represents an operator form that can be applied to a string. Equivalent to StringContainsQ\["string", patt ~~ EndOfString\]. Options: IgnoreCase -\> True treats upper/lowercase as equivalent.

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringEndsQ["commit", "it"]
Out[1]= True

In[2]:= StringEndsQ["commit", "mi"]
Out[2]= False

In[3]:= StringEndsQ[{"apple", "banana"}, "a"]
Out[3]= {False, True}

In[4]:= StringEndsQ["a123", DigitCharacter ..]
Out[4]= True
```

### Applications (3)

A matching suffix

```mathematica
In[5]:= StringEndsQ["report.pdf", ".pdf"]
Out[5]= True
```

No match

```mathematica
In[6]:= StringEndsQ["report.txt", ".pdf"]
Out[6]= False
```

Threads over a list

```mathematica
In[7]:= StringEndsQ[{"a.c", "b.h", "d.c"}, ".c"]
Out[7]= {True, False, True}
```

## Implementation notes

**Algorithm.** `StringEndsQ` is the `SQ_ENDS` face of the shared `sq_dispatch` core — the mirror of `StringStartsQ`. It wraps the pattern as `StringExpression[patt, EndOfString]`, where the translator renders `EndOfString` as `\z` (an absolute end anchor), and hands that unanchored to the shared rule builder before the shared `sq_match` early-exits on the first hit.

**Data structures.** The synthesised wrapper `Expr` (outliving the rule set, as `RegexRule.lhs` borrows into it) and the `RegexRule` array, whole-match pair only.

**Complexity / limits.** Equivalent to `StringContainsQ[s, patt ~~ EndOfString]`. `IgnoreCase` option; curried operator form `StringEndsQ[patt]`. A non-string subject leaves the call unevaluated; byte semantics.

**Attributes:** `Protected`.

## References

- Source: [`src/strings/regex/stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/stringcontainsq.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)

## Notes & additional examples

### Notes

`StringEndsQ[s, patt]` is equivalent to
`StringContainsQ[s, patt ~~ EndOfString]`: it anchors the pattern to the end of
the string by wrapping it as `StringExpression[patt, EndOfString]` before
matching.

It shares the predicate core with `StringContainsQ`, `StringFreeQ`, and
`StringStartsQ`, takes an `IgnoreCase` option, and offers the curried operator
form `StringEndsQ[patt]`.
