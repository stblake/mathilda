/*
 * regex_engine.c - PCRE2 wrapper implementation (see regex_engine.h).
 *
 * The whole PCRE2-facing body lives under #ifdef USE_REGEX; the #else branch
 * provides link-compatible stubs so a build without PCRE2 still succeeds.
 */

#include "regex_engine.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#ifdef USE_REGEX

/* The code-unit width must be fixed before including <pcre2.h>. The makefile /
 * CMake pass -DPCRE2_CODE_UNIT_WIDTH=8; guard here so the header is usable even
 * if a consumer forgets the flag. Mathilda strings are byte-oriented, so the
 * 8-bit library is the correct choice. */
#ifndef PCRE2_CODE_UNIT_WIDTH
#define PCRE2_CODE_UNIT_WIDTH 8
#endif
#include <pcre2.h>

struct RegexProgram {
    pcre2_code*       code;
    pcre2_match_data* md;
    uint32_t          ngroups;   /* capture count, excluding group 0 */
};

int regex_available(void) { return 1; }

RegexProgram* regex_compile(const char* pattern, char* errbuf, size_t errcap) {
    if (!pattern) {
        if (errbuf && errcap) snprintf(errbuf, errcap, "null pattern");
        return NULL;
    }

    int errornumber = 0;
    PCRE2_SIZE erroroffset = 0;
    pcre2_code* code = pcre2_compile((PCRE2_SPTR)pattern,
                                     PCRE2_ZERO_TERMINATED,
                                     0, /* default options */
                                     &errornumber, &erroroffset, NULL);
    if (!code) {
        if (errbuf && errcap) {
            PCRE2_UCHAR msg[256];
            pcre2_get_error_message(errornumber, msg, sizeof(msg));
            snprintf(errbuf, errcap, "%s (at offset %zu)",
                     (const char*)msg, (size_t)erroroffset);
        }
        return NULL;
    }

    /* Best-effort JIT; correctness does not depend on it succeeding. */
    (void)pcre2_jit_compile(code, PCRE2_JIT_COMPLETE);

    pcre2_match_data* md = pcre2_match_data_create_from_pattern(code, NULL);
    if (!md) {
        pcre2_code_free(code);
        if (errbuf && errcap) snprintf(errbuf, errcap, "out of memory");
        return NULL;
    }

    RegexProgram* prog = malloc(sizeof(*prog));
    if (!prog) {
        pcre2_match_data_free(md);
        pcre2_code_free(code);
        if (errbuf && errcap) snprintf(errbuf, errcap, "out of memory");
        return NULL;
    }

    uint32_t cc = 0;
    pcre2_pattern_info(code, PCRE2_INFO_CAPTURECOUNT, &cc);
    prog->code = code;
    prog->md = md;
    prog->ngroups = cc;
    return prog;
}

void regex_free(RegexProgram* prog) {
    if (!prog) return;
    if (prog->md) pcre2_match_data_free(prog->md);
    if (prog->code) pcre2_code_free(prog->code);
    free(prog);
}

int regex_group_count(const RegexProgram* prog) {
    return prog ? (int)prog->ngroups : 0;
}

int regex_match(RegexProgram* prog, const char* subj, size_t len, size_t start,
                size_t* ov, size_t ov_pairs) {
    if (!prog || !subj || !ov || ov_pairs == 0) return -1;
    if (start > len) return 0;

    int rc = pcre2_match(prog->code, (PCRE2_SPTR)subj, (PCRE2_SIZE)len,
                         (PCRE2_SIZE)start, 0, prog->md, NULL);
    if (rc == PCRE2_ERROR_NOMATCH) return 0;
    if (rc < 0) return -1;
    /* rc == 0 means the ovector was too small for all groups; the whole match
     * (pair 0) is still valid, and we only fill up to ov_pairs anyway. */

    PCRE2_SIZE* pov = pcre2_get_ovector_pointer(prog->md);
    size_t avail = (size_t)prog->ngroups + 1;      /* group 0..ngroups */
    for (size_t i = 0; i < ov_pairs; i++) {
        if (i < avail && pov[2 * i] != PCRE2_UNSET) {
            ov[2 * i]     = (size_t)pov[2 * i];
            ov[2 * i + 1] = (size_t)pov[2 * i + 1];
        } else {
            ov[2 * i]     = REGEX_UNSET;
            ov[2 * i + 1] = REGEX_UNSET;
        }
    }
    return 1;
}

#else /* !USE_REGEX -------------------------------------------------------- */

struct RegexProgram { int unused; };

int regex_available(void) { return 0; }

RegexProgram* regex_compile(const char* pattern, char* errbuf, size_t errcap) {
    (void)pattern;
    if (errbuf && errcap)
        snprintf(errbuf, errcap,
                 "regular-expression support not compiled in "
                 "(install PCRE2 and rebuild with USE_REGEX=1)");
    return NULL;
}

void regex_free(RegexProgram* prog) { (void)prog; }
int  regex_group_count(const RegexProgram* prog) { (void)prog; return 0; }

int regex_match(RegexProgram* prog, const char* subj, size_t len, size_t start,
                size_t* ov, size_t ov_pairs) {
    (void)prog; (void)subj; (void)len; (void)start; (void)ov; (void)ov_pairs;
    return -1;
}

#endif /* USE_REGEX */
