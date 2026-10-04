---
source: src/strings/regex/stringreplace.c
---
**Algorithm.** `builtin_stringreplace` builds the rule set *unanchored* and requires every element to carry a replacement RHS (otherwise the call is left unevaluated). `sr_scalar_str` then scans left to right: at each position it takes the leftmost match across all rules (ties broken by rule order), copies the literal text before it, appends the rule's expanded replacement (`regex_rule_replacement`, with `$0`/`$n` substituted), and resumes at the match end. A zero-width match copies one character and advances by one so the loop always makes progress.

**Data structures.** A growable byte buffer (`RegexBuf`); a per-match ovector of up to `REGEX_MAX_PAIRS` pairs holds the whole match and its capture groups for template expansion.

**Complexity / limits.** One left-to-right pass with a PCRE2 probe per rule per position. A list of subjects threads, with non-string elements copied through. Offsets are byte offsets. Only string replacements with `$n` templates are expanded; a non-string RHS falls back to copying the matched text.
