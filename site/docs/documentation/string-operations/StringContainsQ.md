# StringContainsQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringContainsQ["string", patt]`**

Gives True if any substring of "string" matches the string expression patt, and False otherwise.

**`StringContainsQ["string", {p1, p2, ...}]`**

Gives True if any substring matches any of the pi.

**`StringContainsQ[{s1, s2, ...}, patt]`**

Gives the list of results for each of the si.

**`StringContainsQ[patt]`**

Represents an operator form that can be applied to a string. Equivalent to !StringFreeQ\["string", patt\], and to StringMatchQ\["string", \_\_\_ ~~ patt ~~ \_\_\_\]. Options: IgnoreCase -\> True treats upper/lowercase as equivalent.

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringContainsQ["bcde", "b" ~~ __ ~~ "e"]
Out[1]= True

In[2]:= StringContainsQ[{"a", "b", "ab", "abcd", "bcde"}, "a"]
Out[2]= {True, False, True, True, False}

In[3]:= StringContainsQ["bac 123", RegularExpression["a.*"] ~~ DigitCharacter ..]
Out[3]= True

In[4]:= Select[{"abc", "xyz", "bat"}, StringContainsQ["a"]]
Out[4]= {"abc", "bat"}
```

### Options (1)

```mathematica
In[5]:= StringContainsQ["abcd", "BC", IgnoreCase -> True]
Out[5]= True
```

### Worked examples (1)

```mathematica
In[6]:= Options[StringContainsQ]
Out[6]= {IgnoreCase -> False}
```

### Applications (3)

Some substring matches

```mathematica
In[7]:= StringContainsQ["abcdef", "cd"]
Out[7]= True
```

Case folded

```mathematica
In[8]:= StringContainsQ["Mathilda", "math", IgnoreCase -> True]
Out[8]= True
```

Operator form as a predicate

```mathematica
In[9]:= Select[{"cat", "dog", "cow"}, StringContainsQ["o"]]
Out[9]= {"dog", "cow"}
```

## Implementation notes

**Algorithm.** The four substring predicates — `StringContainsQ`, `StringFreeQ`, `StringStartsQ`, `StringEndsQ` — share one core, `sq_dispatch`, and differ only by an `SqKind` tag. It seeds `IgnoreCase` from the registered `Options`, strips trailing option rules, and — with a single positional argument — returns the curried operator form `Function[head[#1, patt, opts...]]`. Otherwise it builds the rule set unanchored and asks `sq_match`: one `regex_match` per rule with an early exit on the first hit (a predicate never needs to know *where* or *how often*). `StringContainsQ` (`SQ_CONTAINS`) reports that match directly.

**Data structures.** A `RegexRule` array; only the whole-match pair (`ov[2]`) is requested — no `regex_scan` enumeration and no capture pool.

**Complexity / limits.** One PCRE2 probe per rule, short-circuited. It is equivalent to `!StringFreeQ[s, p]` and to `StringMatchQ[s, ___ ~~ patt ~~ ___]`. A non-string subject — or a non-string element of a subject list — leaves the *whole* call unevaluated rather than threading a wrong Boolean. Byte semantics.

**Attributes:** `Protected`.

## References

**See also:** [StringMatchQ](../../string-operations/StringMatchQ/), [SetOptions](../../assignment-and-rules/SetOptions/)

- Source: [`src/strings/regex/stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/stringcontainsq.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_parallelmixedspecial.c`](https://github.com/stblake/mathilda/blob/main/tests/test_parallelmixedspecial.c)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)

## Notes & additional examples

### Notes

`StringContainsQ` searches *inside* the string, where `StringMatchQ` anchors to
the whole of it: it is equivalent to `StringMatchQ[s, ___ ~~ patt ~~ ___]` and to
`!StringFreeQ[s, patt]`.

The one-argument operator form `StringContainsQ[patt]` curries into a predicate
(a curried `IgnoreCase -> True` is carried through), which makes it a natural
second argument to `Select`. It shares a single core with `StringFreeQ`,
`StringStartsQ`, and `StringEndsQ`.
