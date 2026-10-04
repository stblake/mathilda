### Worked examples

```mathematica
In[1]:= StringCases["x1y22z333", RegularExpression["\\d+"]]  (* a run of digits *)
```

```mathematica
In[1]:= StringReplace["John Smith", RegularExpression["(\\w+) (\\w+)"] -> "$2 $1"]  (* $n back-references a captured group *)
```

```mathematica
In[1]:= StringMatchQ["abc123", RegularExpression["[a-z]+\\d+"]]  (* letters then digits, whole string *)
```

```mathematica
In[1]:= StringSplit["a1b2c3", RegularExpression["\\d"]]  (* the digits are the delimiters *)
```

### Notes

`RegularExpression["re"]` is an inert head carrying raw PCRE2 source; the string
functions compile it once with `pcre2_compile` and reuse it. In a replacement
right-hand side `$0` is the whole match, `$n` the n-th capture group, and `$$` a
literal `$`.

Backslashes are doubled inside a Mathilda string literal (`\\d` is the single
regex token `\d`). An invalid pattern emits `RegularExpression::regex` but stays
inert; a build without PCRE2 warns and leaves the call unevaluated.
