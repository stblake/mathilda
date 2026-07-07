#ifndef MATHILDA_REGEX_ENGINE_H
#define MATHILDA_REGEX_ENGINE_H

/*
 * regex_engine.h - thin, Expr-agnostic wrapper around PCRE2.
 *
 * This is the ONLY translation unit in the project that includes <pcre2.h>.
 * The wrapper is compiled unconditionally; when Mathilda is built without
 * PCRE2 (USE_REGEX undefined) every entry point degrades to a stub so the
 * rest of the string subsystem links and runs, with regex_available()
 * returning 0 and regex_compile() returning NULL.  This mirrors the
 * USE_MPFR / USE_LAPACK / USE_GRAPHICS graceful-degrade policy.
 *
 * Byte semantics: all offsets are byte offsets into the (raw char*) subject.
 * Mathilda strings are byte-oriented, so this matches the rest of src/strings.
 */

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque compiled program; owns the underlying pcre2_code + match data. */
typedef struct RegexProgram RegexProgram;

/* Returns 1 when built with PCRE2 support, 0 otherwise. */
int regex_available(void);

/*
 * Compile `pattern` (a NUL-terminated PCRE regular expression).
 *
 * On success returns a non-NULL owned program; free it with regex_free().
 * On failure (syntax error, OOM, or built without USE_REGEX) returns NULL and,
 * when errbuf != NULL and errcap > 0, writes a NUL-terminated diagnostic into
 * errbuf.
 */
RegexProgram* regex_compile(const char* pattern, char* errbuf, size_t errcap);

/* Free a program returned by regex_compile(). NULL is ignored. */
void regex_free(RegexProgram* prog);

/* Number of capture groups (excluding whole-match group 0). */
int regex_group_count(const RegexProgram* prog);

/*
 * Attempt one match at or after byte offset `start` in subj[0..len).
 *
 * `ov` receives (start, end) byte offsets for group 0..ngroups, i.e. it must
 * have room for at least ov_pairs*2 size_t entries; unset/absent groups get
 * (REGEX_UNSET, REGEX_UNSET).  Returns 1 on match, 0 on no match, -1 on error.
 */
int regex_match(RegexProgram* prog, const char* subj, size_t len, size_t start,
                size_t* ov, size_t ov_pairs);

/* Sentinel written into `ov` for an unset capture group. */
#define REGEX_UNSET ((size_t)-1)

#ifdef __cplusplus
}
#endif

#endif /* MATHILDA_REGEX_ENGINE_H */
