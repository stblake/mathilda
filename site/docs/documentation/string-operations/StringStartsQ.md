# StringStartsQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringStartsQ["string", patt]`**

Gives True if a prefix of "string" matches the string expression patt, and False otherwise.

**`StringStartsQ["string", {p1, p2, ...}]`**

Gives True if a prefix matches any of the pi.

**`StringStartsQ[{s1, s2, ...}, patt]`**

Gives the list of results for each of the si.

**`StringStartsQ[patt]`**

Represents an operator form that can be applied to a string. Equivalent to StringContainsQ\["string", StartOfString ~~ patt\]. Options: IgnoreCase -\> True treats upper/lowercase as equivalent.

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringStartsQ["commit", "co"]
Out[1]= True

In[2]:= StringStartsQ["commit", "om"]
Out[2]= False

In[3]:= StringStartsQ[{"apple", "banana"}, "a"]
Out[3]= {True, False}

In[4]:= StringStartsQ["a123", LetterCharacter ~~ DigitCharacter ..]
Out[4]= True
```

### Applications (3)

A matching prefix

```mathematica
In[5]:= StringStartsQ["https://x", "http"]
Out[5]= True
```

No match

```mathematica
In[6]:= StringStartsQ["ftp://x", "http"]
Out[6]= False
```

Threads over a list

```mathematica
In[7]:= StringStartsQ[{"abc", "xyz"}, "a"]
Out[7]= {True, False}
```

## Implementation notes

**Algorithm.** `StringStartsQ` is the `SQ_STARTS` face of the shared `sq_dispatch` core. Rather than thread a new anchor mode through the shared rule builder, it synthesises a temporary `StringExpression[StartOfString, patt]` and hands that *unanchored* to `regex_rules_build_ex`; the translator renders `StartOfString` as `\A`, an absolute anchor that pins the match to offset 0 even in an unanchored search. The translator's `group_join` wraps every child in `(?:...)`, so an alternation-bearing pattern anchors correctly as `(?:(?:\A)(?:a|b))` rather than the broken `\Aa|b`. Then the shared `sq_match` early-exits on the first hit.

**Data structures.** The synthesised wrapper `Expr` — which must outlive the rule set, since `RegexRule.lhs` borrows into it — plus the `RegexRule` array (whole-match pair only).

**Complexity / limits.** Equivalent to `StringContainsQ[s, StartOfString ~~ patt]`. `IgnoreCase` option; curried operator form `StringStartsQ[patt]`. A non-string subject leaves the call unevaluated; byte semantics.

**Attributes:** `Protected`.

## References

- Source: [`src/strings/regex/stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/stringcontainsq.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)

## Notes & additional examples

### Notes

`StringStartsQ[s, patt]` is equivalent to
`StringContainsQ[s, StartOfString ~~ patt]`: it anchors the pattern to the start
of the string by wrapping it as `StringExpression[StartOfString, patt]` before
matching.

It shares the predicate core with `StringContainsQ`, `StringFreeQ`, and
`StringEndsQ`, takes an `IgnoreCase` option, and offers the curried operator form
`StringStartsQ[patt]`.
