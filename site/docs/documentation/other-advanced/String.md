# String

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`String`**

is the head of string objects. As a type specification in Read and ReadList it reads a line, up to a newline.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

String is the head of string objects

```mathematica
In[1]:= Head["abc"]
Out[1]= String
```

So _String matches any string

```mathematica
In[2]:= MatchQ["abc", _String]
Out[2]= True
```

The bare symbol String is itself a Symbol

```mathematica
In[3]:= Head[String]
Out[3]= Symbol
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Characters of 200k chars | 4.38 s | 2.54 s | 0.419 s |
| StringSplit on space, 200k chars | 3.12 s | 4.13 s | 0.723 s |
| StringCases regex, 200k chars | 2.36 s | 3.69 s | 3.34 s |
| StringReplace regex, 200k chars | 1.97 s | 4.74 s | 3.35 s |
| StringReplace literal, 200k chars | 0.332 s | 1.17 s | 0.195 s |
| StringCount substring, 200k chars | 0.224 s | 0.364 s | 0.102 s |

## Implementation notes

**Definition.** `String` has two roles. Primarily it is the **head of string
objects**: a string leaf is an `EXPR_STRING`, and `Head` of one returns the symbol
`String`, so `_String` is the pattern that matches any string. Secondarily it is a
type specification in `Read` and `ReadList`, where it reads one line (up to a
newline). The symbol itself is inert — no builtin, no value; its docstring is in
`info.c`.

**Representation.** Strings are stored as `EXPR_STRING` nodes carrying a C string,
not as `String[...]` compounds; `String` is the *head* those leaves report, supplied
by the `Head` builtin in `src/core.c` (the `EXPR_STRING` case), which is why
`StringQ`/`MatchQ[..., _String]` work. As a read type, the scanner in `src/io/read.c`
maps `n == SYM_String` to the `RT_STRING` read tag, whose reader returns the next
line as a string. The bare symbol `String` is itself an `EXPR_SYMBOL`, so
`Head[String]` is `Symbol`.

**Usage & limits.** The head role is pervasive — any pattern, type test, or
`Cases`/`Select` over strings leans on it. The read-type role needs a file or open
stream, so it cannot be exercised on a bare literal. `String` performs no
computation on its own in either role.

**Attributes:** none registered.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`String` is the head carried by string leaves — `Head["text"]` is `String`, which
is what makes `_String` match any string in a pattern, type test, or
`Cases`/`Select`. Strings are stored as `EXPR_STRING` nodes, not as `String[...]`
compounds, so `String` appears as a head but never wraps anything. It is *also* a
`Read`/`ReadList` type specification that reads one line (up to a newline); that
role needs a file or open stream to exercise. The bare symbol `String` is an
ordinary `Symbol`.
