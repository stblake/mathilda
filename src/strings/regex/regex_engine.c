/*
 * regex_engine.c - PCRE2 wrapper implementation (see regex_engine.h).
 *
 * The whole PCRE2-facing body lives under #ifdef USE_REGEX; the #else branch
 * provides link-compatible stubs so a build without PCRE2 still succeeds.
 */

#include "regex_engine.h"
#include "message.h"   /* mth_message: an engine limit */

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
    pcre2_code*          code;
    pcre2_match_data*    md;       /* reused by non-nested matches */
    pcre2_match_context* mctx;     /* NULL unless a callout hook is set */
    RegexCalloutFn       callout;
    void*                callout_data;
    uint32_t             ngroups;  /* capture count, excluding group 0 */
    uint32_t             match_opts;
    int                  busy;     /* a match on this program is in progress */
};

int regex_available(void) { return 1; }

/* PCRE2 callout trampoline: present the block as plain offsets. */
static int callout_trampoline(pcre2_callout_block* cb, void* data) {
    RegexProgram* prog = (RegexProgram*)data;
    size_t ov[REGEX_MAX_PAIRS * 2];
    size_t npairs = (size_t)prog->ngroups + 1;
    if (npairs > REGEX_MAX_PAIRS) npairs = REGEX_MAX_PAIRS;
    ov[0] = (size_t)cb->start_match;
    ov[1] = (size_t)cb->current_position;
    for (size_t i = 1; i < npairs; i++) {
        if (i < cb->capture_top && cb->offset_vector[2 * i] != PCRE2_UNSET) {
            ov[2 * i]     = (size_t)cb->offset_vector[2 * i];
            ov[2 * i + 1] = (size_t)cb->offset_vector[2 * i + 1];
        } else {
            ov[2 * i] = ov[2 * i + 1] = REGEX_UNSET;
        }
    }
    return prog->callout((int)cb->callout_number, (const char*)cb->subject,
                         (size_t)cb->subject_length, ov, npairs,
                         prog->callout_data);
}

RegexProgram* regex_compile_ex(const char* pattern, unsigned flags,
                               RegexCalloutFn callout, void* callout_data,
                               char* errbuf, size_t errcap) {
    if (!pattern) {
        if (errbuf && errcap) snprintf(errbuf, errcap, "null pattern");
        return NULL;
    }

    uint32_t opts = 0, mopts = 0;
    if (flags & REGEX_CASELESS)    opts |= PCRE2_CASELESS;
    if (flags & REGEX_DOTALL)      opts |= PCRE2_DOTALL;
    if (flags & REGEX_UTF)         opts |= PCRE2_UTF | PCRE2_UCP | PCRE2_MATCH_INVALID_UTF;
    if (flags & REGEX_ANCHORED)    opts |= PCRE2_ANCHORED;
    if (flags & REGEX_ENDANCHORED) opts |= PCRE2_ENDANCHORED;
    if (callout)                   opts |= PCRE2_NO_AUTO_POSSESS;

    int errornumber = 0;
    PCRE2_SIZE erroroffset = 0;
    pcre2_code* code = pcre2_compile((PCRE2_SPTR)pattern, PCRE2_ZERO_TERMINATED,
                                     opts, &errornumber, &erroroffset, NULL);
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

    RegexProgram* prog = calloc(1, sizeof(*prog));
    pcre2_match_data* md = pcre2_match_data_create_from_pattern(code, NULL);
    if (!prog || !md) {
        free(prog);
        if (md) pcre2_match_data_free(md);
        pcre2_code_free(code);
        if (errbuf && errcap) snprintf(errbuf, errcap, "out of memory");
        return NULL;
    }
    if (callout) {
        prog->mctx = pcre2_match_context_create(NULL);
        if (!prog->mctx) {
            pcre2_match_data_free(md);
            pcre2_code_free(code);
            free(prog);
            if (errbuf && errcap) snprintf(errbuf, errcap, "out of memory");
            return NULL;
        }
        pcre2_set_callout(prog->mctx, callout_trampoline, prog);
    }

    uint32_t cc = 0;
    pcre2_pattern_info(code, PCRE2_INFO_CAPTURECOUNT, &cc);
    prog->code = code;
    prog->md = md;
    prog->ngroups = cc;
    prog->callout = callout;
    prog->callout_data = callout_data;
    prog->match_opts = mopts;
    return prog;
}

RegexProgram* regex_compile(const char* pattern, char* errbuf, size_t errcap) {
    return regex_compile_ex(pattern, 0, NULL, NULL, errbuf, errcap);
}

void regex_free(RegexProgram* prog) {
    if (!prog) return;
    if (prog->md) pcre2_match_data_free(prog->md);
    if (prog->mctx) pcre2_match_context_free(prog->mctx);
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

    /* A nested match on this program (from inside one of its own callouts)
     * gets private match data. */
    pcre2_match_data* md = prog->md;
    if (prog->busy) {
        md = pcre2_match_data_create_from_pattern(prog->code, NULL);
        if (!md) return -1;
    }
    prog->busy++;
    int rc = pcre2_match(prog->code, (PCRE2_SPTR)subj, (PCRE2_SIZE)len,
                         (PCRE2_SIZE)start, prog->match_opts, md, prog->mctx);
    /* The JIT runs on a fixed 32 kB machine stack, which a backtracking
     * pattern such as (a|b)* exhausts on a subject of a few kB.  The
     * interpreter keeps its backtracking frames on the heap, so rerun there
     * rather than failing the whole call. */
    if (rc == PCRE2_ERROR_JIT_STACKLIMIT)
        rc = pcre2_match(prog->code, (PCRE2_SPTR)subj, (PCRE2_SIZE)len,
                         (PCRE2_SIZE)start, prog->match_opts | PCRE2_NO_JIT, md, prog->mctx);
    prog->busy--;

    int ret;
    if (rc == PCRE2_ERROR_NOMATCH || rc == PCRE2_ERROR_CALLOUT) {
        ret = 0;
    } else if (rc < 0) {
        /* A resource limit is not a "no match": the answer is unknown, and
         * the caller leaves its call unevaluated -- say why, through the
         * funnel (Quiet / Check see it). (a+)+$ against a run of a's ending
         * in b backtracks exponentially. */
        if (rc == PCRE2_ERROR_MATCHLIMIT || rc == PCRE2_ERROR_DEPTHLIMIT ||
            rc == PCRE2_ERROR_HEAPLIMIT)
            mth_message("RegularExpression", "limit",
                        "The pattern matcher's backtracking limit was reached; "
                        "the match could not be decided.");
        ret = -1;
    } else {
        /* rc == 0 means the ovector was too small for all groups; the whole
         * match (pair 0) is still valid, and we only fill up to ov_pairs. */
        PCRE2_SIZE* pov = pcre2_get_ovector_pointer(md);
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
        ret = 1;
    }
    if (md != prog->md) pcre2_match_data_free(md);
    return ret;
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

RegexProgram* regex_compile_ex(const char* pattern, unsigned flags,
                               RegexCalloutFn callout, void* callout_data,
                               char* errbuf, size_t errcap) {
    (void)flags; (void)callout; (void)callout_data;
    return regex_compile(pattern, errbuf, errcap);
}

void regex_free(RegexProgram* prog) { (void)prog; }
int  regex_group_count(const RegexProgram* prog) { (void)prog; return 0; }

int regex_match(RegexProgram* prog, const char* subj, size_t len, size_t start,
                size_t* ov, size_t ov_pairs) {
    (void)prog; (void)subj; (void)len; (void)start; (void)ov; (void)ov_pairs;
    return -1;
}

#endif /* USE_REGEX */
