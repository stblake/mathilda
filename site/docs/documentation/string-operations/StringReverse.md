# StringReverse

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringReverse["string"]`**

Reverses the order of the characters in "string".

**`StringReverse[{s1, s2, ...}]`**

Gives the list of results for each of the si. StringReverse is Listable, so it threads automatically over lists.

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= StringReverse["abcdef"]
Out[1]= "fedcba"

In[2]:= StringReverse[{"cat", "dog", "fish", "coelenterate"}]
Out[2]= {"tac", "god", "hsif", "etaretneleoc"}

In[3]:= StringReverse[""]
Out[3]= ""

In[4]:= StringReverse[x]
Out[4]= StringReverse[x]

In[5]:= StringReverse[] StringReverse::argx: StringReverse called with 0 arguments; 1 argument is expected.
```

### Applications (3)

Reverse the characters

```mathematica
In[6]:= StringReverse["abcdef"]
Out[6]= "fedcba"
```

A palindrome is its own reverse

```mathematica
In[7]:= StringReverse["racecar"]
Out[7]= "racecar"
```

Listable, so a list threads

```mathematica
In[8]:= StringReverse[{"ab", "cd"}]
Out[8]= {"ba", "dc"}
```

## Implementation notes

**Algorithm.** `builtin_stringreverse` requires one argument (else `StringReverse::argx`) and copies the subject's bytes into a fresh buffer in reverse order. `StringReverse` is `Listable`, so the evaluator threads it element-wise over a list before the builtin runs; each call therefore sees a single string.

**Data structures.** One output buffer of the same length as the input.

**Complexity / limits.** `O(len)`, byte-oriented (characters are single `char`s, no UTF-8 decoding). A single non-string argument leaves the call unevaluated so symbolic arguments flow through unchanged.

**Attributes:** `Listable`, `Protected`.

## References

- Source: [`src/strings/stringreverse.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringreverse.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_strings.c`](https://github.com/stblake/mathilda/blob/main/tests/test_strings.c)

## Notes & additional examples

### Notes

`StringReverse` reverses the bytes of a string into a fresh buffer. It is
`Listable`, so the evaluator threads it over a list before the builtin runs and
each call sees a single string.

A call with a number of arguments other than one emits `StringReverse::argx`; a
non-string argument leaves the call unevaluated.
