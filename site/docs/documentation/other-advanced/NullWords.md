# NullWords

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`NullWords`**

is an option for Read and ReadList specifying whether a null word is assumed between two adjacent word separators. The default is False.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

A bare Read/ReadList option symbol

```mathematica
In[1]:= Head[NullWords]
Out[1]= Symbol
```

The default is False

```mathematica
In[2]:= Options[ReadList, NullWords]
Out[2]= {NullWords -> False}
```

```mathematica
In[3]:= str = OpenWrite["/tmp/mathilda_nw.txt"]; WriteString[str, "a  b\n"]; Close[str];
```

Default: the empty field is skipped

```mathematica
In[4]:= ReadList["/tmp/mathilda_nw.txt", Word]
Out[4]= {"a", "b"}
```

Now kept as ""

```mathematica
In[5]:= ReadList["/tmp/mathilda_nw.txt", Word, NullWords -> True]
Out[5]= {"a", "", "b"}
```

## Implementation notes

**Definition.** `NullWords` is an **option name** for `Read` and `ReadList` specifying
whether a null (empty) word is assumed between two adjacent word separators. It is a bare
option symbol, not a function — no builtin, no DownValues, not even the `Protected`
attribute — carrying only the docstring in `info_init` (`src/info.c`). The reader
`src/io/read.c` reads it out of the option list (its default `NullWords -> False` is part
of `Options[Read]` / `Options[ReadList]`).

**Representation.** `NullWords` stays an inert `EXPR_SYMBOL`, appearing only on the left of
a rule; `FullForm[NullWords -> True]` is `Rule[NullWords, True]`. The value is a boolean;
the default is `False`.

**Usage & limits.** With `NullWords -> False` (the default) the word reader skips an empty
field between adjacent separators — two spaces read as one gap. With `NullWords -> True` it
emits an empty string `""` for that field, so `"a  b"` reads as `{"a", "", "b"}`. It is one
of the tokenisation options alongside `WordSeparators`, `TokenWords`, `RecordSeparators`
and `NullRecords`, all defaulting from `Options[Read]` so `SetOptions[Read, ...]` adjusts
them globally.

**Attributes:** none registered.

## References

- Source: [`src/io/read.c`](https://github.com/stblake/mathilda/blob/main/src/io/read.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`NullWords` is an option for `Read` and `ReadList` deciding whether an empty word is
assumed between two adjacent word separators. It is an inert option-name symbol — no
builtin, no DownValues — read out of a call's option list by the reader.

The default `False` makes the word reader skip the empty field, so the double space above
reads as `{"a", "b"}`; `NullWords -> True` keeps it, giving `{"a", "", "b"}`. It is one of
the tokenisation options (with `WordSeparators`, `TokenWords`, `RecordSeparators`,
`NullRecords`), all defaulting from `Options[Read]`, so `SetOptions[Read, NullWords ->
True]` would change the global default.
