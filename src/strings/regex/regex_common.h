#ifndef MATHILDA_REGEX_COMMON_H
#define MATHILDA_REGEX_COMMON_H

/*
 * regex_common.h - Expr-aware helpers shared by the regex string builtins
 * (StringMatchQ / StringCases / StringReplace / StringSplit).
 *
 * This layer sits above the pure PCRE2 wrapper (regex_engine.h): it turns a
 * Mathilda pattern argument (RegularExpression[...], a literal string, a
 * Rule/RuleDelayed, or a List of those) into an array of compiled rules, and
 * provides $n replacement-template expansion and a small growable string
 * buffer.  It compiles unconditionally; when PCRE2 is absent every build fails
 * cleanly via regex_available().
 */

#include <stddef.h>
#include "expr.h"
#include "regex_engine.h"

/* Largest capture group we ever expose ($0..$63). */
#define REGEX_MAX_PAIRS 64

/*
 * One compiled pattern, optionally paired with a replacement RHS taken from a
 * Rule/RuleDelayed.  `rhs` is BORROWED from the caller's expression (never
 * owned/freed here); it is NULL for a bare pattern.
 */
typedef struct {
    RegexProgram* prog;   /* owned; free via regex_rules_free */
    Expr*         rhs;    /* borrowed; NULL when the element was a bare pattern */
    Expr*         lhs;    /* borrowed; the delimiter pattern (rule LHS or bare pattern).
                             Lets StringSplit bind a named Pattern[x,...] for `:>` RHS. */
} RegexRule;

/*
 * Build an array of RegexRule from a pattern argument.  Accepts
 * RegularExpression["re"], a literal string (matched literally), Rule/
 * RuleDelayed (rhs recorded), or a List of any of those (alternatives / rule
 * set).  When `anchored` is nonzero each pattern is wrapped \A(?:...)\z for a
 * whole-string match (StringMatchQ).
 *
 * Returns the rule count (>= 1) with a malloc'd array stored in *out, or -1 on
 * an unsupported pattern, a compile error, or when regex support is absent
 * (a diagnostic keyed on `head` is printed to stderr in the latter two cases).
 * A successful result must be released with regex_rules_free().
 */
int  regex_rules_build(Expr* patt, int anchored, RegexRule** out, const char* head);
/*
 * Extended form: when `caseless` is nonzero each pattern matches case-insensitively
 * (a `(?i)` inline modifier is prepended to the compiled source). `regex_rules_build`
 * is the thin `caseless == 0` wrapper. Used by StringSplit's IgnoreCase option.
 */
int  regex_rules_build_ex(Expr* patt, int anchored, int caseless,
                          RegexRule** out, const char* head);
void regex_rules_free(RegexRule* rules, int n);

/*
 * Translate a Wolfram string pattern into malloc'd PCRE source (caller frees), or
 * return NULL if the pattern is unsupported. Handles literal strings,
 * RegularExpression["re"], the character-class heads (Whitespace, WhitespaceCharacter,
 * LetterCharacter, DigitCharacter, WordCharacter, NumberString), StringExpression,
 * Alternatives, Repeated/RepeatedNull, Except, Blank, Pattern (capture group), and
 * PatternTest with a known predicate (LetterQ/DigitQ/UpperCaseQ/LowerCaseQ).
 *
 * When `is_null` is non-NULL it is set to 1 iff `patt` is the empty string "" (the
 * "split at every character" null delimiter); callers that care handle it specially.
 */
char* wl_pattern_to_regex(Expr* patt, int* is_null);

/*
 * Expand a replacement template (`$0`..`$N`, and `$$` -> literal `$`) using the
 * capture offsets `ov` (npairs pairs of byte offsets) against subject `subj`.
 * Unknown/unset groups expand to the empty string.  Returns a malloc'd
 * NUL-terminated string (caller frees), or NULL on OOM.
 */
char* regex_expand_template(const char* tpl, const char* subj,
                            const size_t* ov, size_t npairs);

/*
 * Resolve rule `r`'s replacement for one match.  Handles the supported case
 * (rhs is a string, with $n expansion) and returns NULL when the RHS is
 * unsupported in the current scope (a non-string RHS).  Caller frees.
 */
char* regex_rule_replacement(const RegexRule* r, const char* subj,
                             const size_t* ov, size_t npairs);

/* Small growable byte buffer used to assemble result strings. Zero-initialise
 * (RegexBuf b = {0};) and release with free(b.p). */
typedef struct { char* p; size_t len, cap; } RegexBuf;

/* Append n bytes; keeps the buffer NUL-terminated. Returns 0, or -1 on OOM. */
int regexbuf_add(RegexBuf* b, const char* s, size_t n);

#endif /* MATHILDA_REGEX_COMMON_H */
