### Worked examples

```mathematica
In[1]:= Names["Sin"]
```

```mathematica
In[1]:= Names["Ar@"]  (* @ matches a run of one-or-more non-uppercase letters *)
```

```mathematica
In[1]:= Names[RegularExpression["Si."]]  (* a PCRE pattern, anchored to the whole name *)
```

```mathematica
In[1]:= MemberQ[Names["System`*"], "System`Sin"]  (* a backtick pattern enumerates builtins *)
```

### Notes

`Names["patt"]` gives the canonically sorted list of symbol names matching a
string pattern. The pattern is matched against the **whole** name (anchored) with
two metacharacters: `*` for zero or more of anything, and `@` for one or more
characters that are not uppercase letters (so `"Ar@"` stops before the capital in
`ArcSin`). Every other character, the context `` ` `` included, is literal; a
pattern element may instead be `RegularExpression["re"]`.

A pattern containing a backtick is matched against, and returns, each symbol's
fully context-qualified name (`System`Sin`, `Global`x`), which is what makes
`Names["System`*"]` enumerate all builtins. The result is always sorted, so
`Names[patt]` equals `Sort[Names[patt]]`. On a build without PCRE2 a
`RegularExpression` pattern emits `Names::regavail` and stays unevaluated.
