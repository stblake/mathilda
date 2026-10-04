# WordSeparators

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`WordSeparators`**

is an option for Read and ReadList giving the list of strings that separate words. The default is {" ", "\t"}.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

An inert option keyword: it just names the left of a rule

```mathematica
In[1]:= WordSeparators -> {" ", ","}
Out[1]= WordSeparators -> {" ", ","}
```

Its default is {" ", "\t"}, alongside the other read options

```mathematica
In[2]:= Options[ReadList]
Out[2]= {RecordSeparators -> {"\r\n", "\n", "\r"}, WordSeparators -> {" ", "\t"}, TokenWords -> {}, NullRecords -> False, NullWords -> False}
```

A plain inert symbol, with no value of its own

```mathematica
In[3]:= Attributes[WordSeparators]
Out[3]= {}
```

## Implementation notes

**Definition.** `WordSeparators` is an option for `Read` and `ReadList` giving the
list of strings that separate words when the type specification is `Word`. The
default is `{" ", "\t"}`. It is an inert option keyword — no builtin and no value of
its own (it carries no attributes) — declared with its docstring in `info.c` and it
appears in `Options[Read]` / `Options[ReadList]` with its default.

**Representation.** A bare `EXPR_SYMBOL`. During a read, the options are gathered
into a small table (`src/io/read.c`, mirrored in `readlist.c`): the `"WordSeparators"`
row binds a local pointer `o_word` from the user's rules, which `read_cfg_build`
folds into the read configuration. The tokeniser then uses that separator list when
cutting `Word` tokens out of the stream (`RT_WORD`), falling back to the string-list
default when the option is not given.

**Usage & limits.** Meaningful only in a `Read`/`ReadList` call that reads `Word`
tokens; it does not affect `Number`, `Record` (which uses `RecordSeparators`) or
`String` (whole line) reads. The symbol performs no computation — it names the
option whose value the reader consults.

**Attributes:** none registered.

## References

- Source: [`src/io/read.c`](https://github.com/stblake/mathilda/blob/main/src/io/read.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`WordSeparators` is an option for `Read` and `ReadList` giving the strings that
separate words when the type specification is `Word`; the default is `{" ", "\t"}`.
It is an inert keyword — the reader consults its value when cutting `Word` tokens
from a stream. A read that uses it, such as `ReadList["data.txt", Word,
WordSeparators -> {",", " "}]`, needs a file (or open stream) on disk, so the
runnable examples above show it only in isolation — as the left of an option rule
and as the default surfaced by `Options[ReadList]`. It does not affect `Number`,
`Record` (which uses `RecordSeparators`), or whole-line `String` reads.
