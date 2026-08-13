# KeyValuePattern

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyValuePattern[{k1 -> p1, ...}]`**

A pattern matching an association (or list of rules) that contains keys matching k1, ... with values matching p1, .... Value patterns may bind (e.g. KeyValuePattern\[{"a" -\> v\_}\]). KeyValuePattern\[k -\> p\] is the single-key form.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= MatchQ[<|"a" -> 1, "b" -> 2|>, KeyValuePattern[{"a" -> _}]]
Out[1]= True

In[2]:= Replace[<|"a" -> 5, "b" -> 2|>, KeyValuePattern[{"a" -> v_}] :> v]
Out[2]= 5

In[3]:= Cases[{<|"t" -> 1|>, <|"t" -> 2|>, <|"x" -> 3|>}, KeyValuePattern[{"t" -> _}]]
Out[3]= {<|"t" -> 1|>, <|"t" -> 2|>}
```

### Scope (2)

```mathematica
In[4]:= Cases[{<|"p" -> 3|>, <|"p" -> 9|>}, KeyValuePattern[{"p" -> v_}] /; v > 5 :> v]
Out[4]= {9}

In[5]:= area[KeyValuePattern[{"w" -> w_, "h" -> h_}]] := w h; area[<|"w" -> 3, "h" -> 4|>]
Out[5]= 12
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/patterns.c`](https://github.com/stblake/mathilda/blob/main/src/patterns.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
