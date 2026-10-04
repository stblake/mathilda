# KeyValuePattern

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyValuePattern[{k1 -> p1, ...}]`**

A pattern matching an association (or list of rules) that contains keys matching k1, ... with values matching p1, .... Value patterns may bind (e.g. KeyValuePattern\[{"a" -\> v\_}\]). KeyValuePattern\[k -\> p\] is the single-key form.

## Examples (7)

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

### Applications (2)

Has key a with any value

```mathematica
In[6]:= MatchQ[<|a -> 1, b -> 2|>, KeyValuePattern[{a -> _}]]
Out[6]= True
```

Pull the y-values

```mathematica
In[7]:= Cases[{<|x -> 1, y -> 2|>, <|x -> 5|>}, KeyValuePattern[{y -> v_}] :> v]
Out[7]= {2}
```

## Implementation notes

**Algorithm.** `KeyValuePattern` is an inert pattern head consumed by the matcher,
not an evaluation builtin (it is interned and marked `Protected` in
`src/patterns.c`; the matching lives in `match_internal`). `KeyValuePattern[{k1 -> p1,
…}]` matches an association — or a list of rules — that contains, for each
requirement `ki -> pi`, an entry whose key matches `ki` and whose value matches `pi`.
The requirements are solved with backtracking (`kvp_match_reqs`), so a consistent
assignment is found whenever one exists even when requirements share bound variables
(e.g. `KeyValuePattern[{k_ -> _, _ -> k_}]`). An empty spec matches any association,
and a bare single rule is the one-requirement form.

**Data structures.** The requirement list is read directly off the pattern's
argument; bindings accumulate in the matcher's `MatchEnv`, with `env_rollback` undoing
a failed branch of the backtracking search.

**Complexity / limits.** Each requirement is matched against the subject's entries;
the backtracking is exponential only in the number of requirements that share
variables, which is small in practice. The subject must have head `Association` or
`List`, or the match fails immediately.

**Attributes:** `Protected`.

## References

- Source: [`src/match.c`](https://github.com/stblake/mathilda/blob/main/src/match.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_atomicity.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_atomicity.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

`KeyValuePattern` matches an association (or list of rules) that *contains* the given
`key -> pattern` entries, in any order and ignoring extra keys. Value patterns may
bind, so `KeyValuePattern[{"a" -> v_}]` captures the value at `"a"`. Requirements that
share a bound variable are resolved consistently by backtracking.
