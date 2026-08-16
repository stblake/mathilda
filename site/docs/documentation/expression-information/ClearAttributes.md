# ClearAttributes

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ClearAttributes[s, attr] removes attr from the list of attributes of s.`**

**`ClearAttributes[s, {attr1, attr2, ...}] removes several attributes at a time.`**

**`ClearAttributes[{s1, s2, ...}, attrs] removes attributes from several symbols at a time.`**

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (12)

```mathematica
In[1]:= SetAttributes[f, Listable]

In[2]:= f[{1, 2, 3}]
Out[2]= {f[1], f[2], f[3]}

In[3]:= ClearAttributes[f, Listable]

In[4]:= f[{1, 2, 3}]
Out[4]= f[{1, 2, 3}]

In[5]:= SetAttributes[f, {Flat, Orderless, OneIdentity}]

In[6]:= ClearAttributes[f, OneIdentity]

In[7]:= Attributes[f]
Out[7]= {Flat, Orderless}

In[8]:= ClearAttributes[f, {Flat, Orderless}]

In[9]:= Attributes[f]
Out[9]= {}

In[10]:= SetAttributes[{g, h}, Protected]

In[11]:= ClearAttributes[{g, h}, Protected]

In[12]:= Attributes[g]
Out[12]= {}
```

## Implementation notes

`builtin_clear_attributes` (`src/attr.c`) clears the bitflags named in its second argument from the target symbol(s) via `clear_attributes_for_symbol`. The first argument may be one symbol/string or a `List` of them; it returns `Null`. `ClearAttributes` carries `ATTR_HOLDFIRST` so the symbol is not evaluated first.

- `HoldFirst`, `Protected`.
- `ClearAttributes` modifies `Attributes[s]`.
- Cannot clear attributes of a `Locked` symbol.
- Clearing an attribute that is not set is a no-op.

**Attributes:** `HoldFirst`, `Protected`.

## References

**See also:** [HoldFirst](../../other-advanced/HoldFirst/)

- Source: [`src/attr.c`](https://github.com/stblake/mathilda/blob/main/src/attr.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_eval_timestamps.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eval_timestamps.c)
