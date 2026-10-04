---
references:
  - "P. Hazel, *PCRE2 — Perl-Compatible Regular Expressions*, 10.x, https://www.pcre.org/."
source: src/strings/regex/regularexpression.c
---
**Algorithm.** `RegularExpression["re"]` is an inert data head: `builtin_regularexpression` validates a single string argument and otherwise returns `NULL`, so the head survives evaluation unchanged and carries raw PCRE source for the consuming builtins (`StringMatchQ`, `StringCases`, `StringReplace`, `StringSplit`, `StringPosition`, …). As a best-effort syntax check it runs `regex_compile` once when PCRE2 is present, emitting `RegularExpression::regex` on a bad pattern while staying inert.

The real work is in the two shared layers the consumers call. The translator `wl_pattern_to_regex` (`src/strings/regex/string_pattern.c`) passes `RegularExpression["re"]` through verbatim as PCRE source (and renders the symbolic string-pattern heads — `Whitespace`, `LetterCharacter`, `~~`, `|`, `..`, `Except`, `Pattern` → capture/backreference — into PCRE too). The engine wrapper `regex_compile` (`src/strings/regex/regex_engine.c`) calls `pcre2_compile` plus a best-effort `pcre2_jit_compile`.

**Data structures.** A `RegexProgram` wraps a `pcre2_code`, a `pcre2_match_data`, and the capture count; the 8-bit PCRE2 library is used because Mathilda strings are byte-oriented. In a replacement RHS, `$0` is the whole match and `$n` the n-th group, expanded by `regex_expand_template`.

**Complexity / limits.** Matching cost is PCRE2's; up to `$0..$63` capture pairs are exposed (`REGEX_MAX_PAIRS`). Backslashes are doubled in a Mathilda string literal (`\\d` → `\d`). Without PCRE2 (`USE_REGEX` unset) the engine stubs out and every string-pattern builtin warns (`regavail`) and stays unevaluated.
