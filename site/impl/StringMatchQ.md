---
source: src/strings/regex/stringmatchq.c
---
**Algorithm.** `builtin_stringmatchq` compiles the pattern *anchored* — `regex_rules_build` with `anchored = 1` wraps each rule as `\A(?:...)\z` — so a match means the pattern covers the whole subject. For each subject it runs `regex_match` at offset 0 and returns `True` on the first rule that matches, else `False`. A `List` pattern becomes several rules (a match if any one matches); a list of subjects threads, with non-string elements copied through unchanged.

The pattern may be a literal string, `RegularExpression[...]`, a general string expression (`~~`, `|`, `..`, character classes, `Pattern`), or a list of alternatives — all handled by the shared translator `wl_pattern_to_regex`.

**Data structures.** A `RegexRule` array, each a compiled `RegexProgram`. Only the whole-match pair is requested (`ov[2]`): a yes/no answer needs no capture pool and no `regex_scan` enumeration.

**Complexity / limits.** One PCRE2 match per rule per subject, early-exiting on the first hit. Byte semantics throughout (no UTF-8 decoding). A non-string subject leaves the call unevaluated.
