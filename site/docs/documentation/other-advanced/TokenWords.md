# TokenWords

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`TokenWords`**

is an option for Read and ReadList giving strings to be read as separate words even when not surrounded by word separators. The default is {}.

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (1)

"+" read as its own word

```mathematica
In[1]:= s = OpenWrite["/tmp/mathilda_tok.txt"]; WriteString[s, "a+b+c"]; Close[s]; ReadList["/tmp/mathilda_tok.txt", Word, TokenWords -> {"+"}]
Out[1]= {"a", "+", "b", "+", "c"}
```

## Implementation notes

**Definition.** `TokenWords` is an **option** for `Read` and `ReadList` giving a list of
strings that are to be read as separate words even when they are not surrounded by word
separators. It has no builtin and no value of its own — an inert, `Protected` keyword
carried in the reader's option list, with default `{}`.

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_TokenWords`). Both `read.c` and
`readlist.c` seed their default option set with `rule_of(SYM_TokenWords, {})`. The
reader's `next_word` scanner checks, at each position, whether a `TokenWords` string
begins there (`token_len`); if so it emits that token as its own field and advances past
it, so an operator glued to its operands (`a+b` with `TokenWords -> {"+"}`) splits into
`a`, `+`, `b`.

**Usage & limits.** `Protected`; meaningful only as a `Read`/`ReadList` option for the
`Word` read type. It has a registered default, so `SetOptions[Read, TokenWords -> ...]`
works. Related tokenisation options are `WordSeparators`, `RecordSeparators`,
`NullWords` and `NullRecords`.

**Attributes:** none registered.

## References

- Source: [`src/io/read.c`](https://github.com/stblake/mathilda/blob/main/src/io/read.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`TokenWords` is an option for `Read` and `ReadList` giving strings that are read as
separate words even when they are not surrounded by word separators. The default is `{}`.

With `TokenWords -> {"+"}`, the `Word` reader emits `+` as its own field, so `"a+b+c"`
reads as `{"a", "+", "b", "+", "c"}` — the standard way to tokenise an expression whose
operators are glued to their operands. It is an inert, `Protected` keyword with a
registered default, so `SetOptions[Read, TokenWords -> ...]` works. Related tokenisation
options are `WordSeparators`, `RecordSeparators`, `NullWords` and `NullRecords`.
