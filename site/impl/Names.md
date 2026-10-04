---
source: src/names.c
---
**Algorithm.** `builtin_names` (`src/names.c`) enumerates the symbol table by
pattern. Each argument element is compiled by `build_pat` into either a glob
(`EXPR_STRING`) or a PCRE2 program (`RegularExpression["re"]`, anchored as
`\A(?:re)\z`); a `List` argument becomes a set of alternative patterns, and no
argument means "match every name". The glob matcher `wl_glob_match` is a recursive
backtracking whole-string match with two metacharacters — `*` (zero or more of
anything) and `@` (one or more non-uppercase characters) — every other byte,
including the context `` ` ``, literal. `symtab_for_each` then visits every symbol,
`match_one` decides the emitted string, and the collected names are `qsort`ed by
`expr_compare` so that `Names[p]` is identical to `Sort[Names[p]]`.

**Data structures.** A `NameCollect` growable `char**` buffer accumulates owned
name strings during the visitor callback; each becomes an `EXPR_STRING` leaf of the
returned `List`. Context qualification is lazy: `full_qualified_name` prefixes a
bare name with `System`` when its `SymbolDef` has a `builtin_func` (or is a
kernel-interned System symbol) and `Global`` otherwise, caching the result across
the pattern loop for one symbol.

**Complexity / limits.** `O(S · P · L)` for `S` symbols, `P` patterns, and glob
cost `L` per name (regex cost dominates when a `RegularExpression` is used). A
backtick in a glob switches both the match target and the returned string to the
fully-qualified name, which is what lets `Names["System`*"]` list the builtins. A
`RegularExpression` pattern on a build without PCRE2 emits `Names::regavail`
(through `mth_message`) and leaves the call unevaluated. Attributes `Protected`.
