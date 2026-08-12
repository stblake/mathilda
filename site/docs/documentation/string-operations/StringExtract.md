# StringExtract

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringExtract["string", n]`**

Extracts the nth whitespace-delimited block of "string" (-n counts from the end).

**`StringExtract["string", spec]`**

Selects blocks with spec: n, -n, {n1, n2, ...}, n1;;n2, or All.

**`StringExtract["string", sep -> pos]`**

Delimits blocks with the string pattern sep. sep -\> All equals StringSplit\["string", sep\].

**`StringExtract["string", pos1, pos2, ...]`**

Extracts across levels: whitespace at the lowest level, then "\n", then "\n\n", and so on for higher levels.

**`StringExtract["string", sep1 -> pos1, sep2 -> pos2, ...]`**

Uses sepi as the separator for successive levels. Absent blocks yield Missing\["PartAbsent", pos\]. A list of strings threads.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringExtract["a bbb  cccc aa   d", 2]
Out[1]= "bbb"

In[2]:= StringExtract["a bbb  cccc aa   d", 2 ;; 4]
Out[2]= {"bbb", "cccc", "aa"}

In[3]:= StringExtract["a--bbb--ccc--dddd", "--" -> 3]
Out[3]= "ccc"

In[4]:= StringExtract["a 1\nb 2\nc 3 x", All, 3]
Out[4]= {Missing["PartAbsent", 3], Missing["PartAbsent", 3], "x"}
```

## Algorithm

stringextract.c - StringExtract[subject, spec_1, ..., spec_k]

Splits a string into blocks and selects blocks by position. Each spec_j is a level: either a bare position (n, -n, {n1,...}, n1;;n2, All) or a rule

```text
`sep -> pos`. Level j splits its input on separator sep_j, then selects with
```

pos_j; if more levels remain it recurses into each selected block.

The split is delegated to StringSplit (built here as StringSplit[str, sep] and evaluated), so the entire WL string-pattern engine, whitespace-run collapsing, and empty-end trimming are reused verbatim. This makes the documented equivalences exact:

```text
    StringExtract[s, patt -> All]  ==  StringSplit[s, patt]
    StringExtract[s, {p1, p2, ...}] == Part[StringSplit[s], {p1, p2, ...}]
```

Bare positions get depth-default separators: the last (lowest) level splits on whitespace, the level above on a single "\n", the next on "\n\n", and so on. An out-of-range single index yields Missing["PartAbsent", n] (matching Wolfram and the Missing["KeyAbsent", ...] convention in part.c), not the usual NULL-unevaluated result.

Registered by regex_init() (regex_init.c) alongside its StringSplit dependency. Docstring lives in info.c. Attribute: Protected.

## Implementation notes

**Attributes:** `Protected`.

## See also

[StringSplit](../../string-operations/StringSplit/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringfns.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringfns.c)
