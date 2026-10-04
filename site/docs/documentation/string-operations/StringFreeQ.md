# StringFreeQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringFreeQ["string", patt]`**

Gives True if no substring of "string" matches the string expression patt, and False otherwise.

**`StringFreeQ["string", {p1, p2, ...}]`**

Gives True if no substring matches any of the pi.

**`StringFreeQ[{s1, s2, ...}, patt]`**

Gives the list of results for each of the si.

**`StringFreeQ[patt]`**

Represents an operator form that can be applied to a string. Equivalent to !StringContainsQ\["string", patt\]. Options: IgnoreCase -\> True treats upper/lowercase as equivalent.

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= StringFreeQ["abcd", "a"]
Out[1]= False

In[2]:= StringFreeQ["abcade", x_ ~~ x_]
Out[2]= True

In[3]:= StringFreeQ[{"ability", "listable", "argument"}, "a" ~~ __ ~~ "t" ~~ ___]
Out[3]= {False, True, False}
```

### Options (1)

```mathematica
In[4]:= StringFreeQ["ac", IgnoreCase -> True]["BACCD"]
Out[4]= False
```

### Applications (3)

Nothing matches, so free

```mathematica
In[5]:= StringFreeQ["abcdef", "xyz"]
Out[5]= True
```

A substring matches, so not free

```mathematica
In[6]:= StringFreeQ["abcdef", "cd"]
Out[6]= False
```

Letters are present

```mathematica
In[7]:= StringFreeQ["hello", LetterCharacter]
Out[7]= False
```

## Implementation notes

**Algorithm.** `StringFreeQ` is the `SQ_FREE` face of the shared `sq_dispatch` core (`src/strings/regex/stringcontainsq.c`). The matching is identical to `StringContainsQ` — one `regex_match` per rule with an early exit — but `sq_answer` inverts the result, so `StringFreeQ` is `True` exactly when no rule matches anywhere. Option seeding, the curried operator form, and list-threading are all shared.

**Data structures.** As `StringContainsQ`: a `RegexRule` array, whole-match pair only, no capture pool.

**Complexity / limits.** `StringFreeQ[s, p] == !StringContainsQ[s, p]` by construction. An `IgnoreCase` option is honoured; a list of patterns means "free of any of them" (an empty list matches nothing). A non-string subject leaves the call unevaluated; byte semantics.

**Attributes:** `Protected`.

## References

**See also:** [StringContainsQ](../../string-operations/StringContainsQ/)

- Source: [`src/strings/regex/stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/stringcontainsq.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)

## Notes & additional examples

### Notes

`StringFreeQ` is the exact negation of `StringContainsQ` and shares its
implementation (one `SqKind` tag selects the inverting answer). It is `True` when
no substring matches the pattern.

It takes an `IgnoreCase` option and offers the curried operator form
`StringFreeQ[patt]`. A list of patterns means "free of any of them".
