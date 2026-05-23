#ifndef MATHILDA_SIMP_INTERNAL_H
#define MATHILDA_SIMP_INTERNAL_H

/*
 * simp_internal.h -- cross-module surface for the simp.c split.
 *
 * Every function declared here is defined in one of the simp_*.c sibling
 * translation units and is referenced from at least one other simp_*.c
 * unit. Functions used only within their own file remain file-static and
 * do NOT appear here. See simp.c for the file-by-file map.
 */

#include "expr.h"
#include "simp.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <time.h>

/* ------------------------------------------------------------------ */
/* Shared types                                                       */
/* ------------------------------------------------------------------ */

/* Candidate set used by the heuristic search loop. Lives in simp_util.c
 * (cs_* helpers) and is consumed by simp_search.c. */
#define SIMP_CAND_CAP 12
#define SIMP_ROUNDS   2

typedef struct {
    Expr** items;
    size_t count;
    size_t capacity;
} CandSet;

/* Per-Simplify memoisation table for the bottom-up driver. Defined and
 * managed in simp_bottomup.c; instantiated as a local in builtin_simplify
 * (simp_builtins.c) and passed by pointer through the bottom-up call
 * chain. */
#define SIMP_MEMO_BUCKETS 256

typedef struct SimpMemoEntry {
    Expr* key;
    Expr* value;
    struct SimpMemoEntry* next;
} SimpMemoEntry;

typedef struct {
    SimpMemoEntry* buckets[SIMP_MEMO_BUCKETS];
} SimpMemo;

/* ------------------------------------------------------------------ */
/* simp_util.c                                                        */
/* ------------------------------------------------------------------ */

Expr* call_unary_owned(const char* head_name, Expr* arg);
Expr* call_unary_copy(const char* head_name, const Expr* arg);
Expr* traced_call_unary(const char* xform, const Expr* in);

bool   simp_debug_enabled(void);
double simp_debug_elapsed_ms(clock_t t0);
void   simp_debug_log(const char* xform, const Expr* in,
                      const Expr* out, double elapsed_ms);

void cs_init(CandSet* cs);
void cs_free(CandSet* cs);
bool cs_contains(const CandSet* cs, const Expr* e);
void cs_add_or_free(CandSet* cs, Expr* e);

size_t score_with_func(const Expr* e, const Expr* complexity_func);

bool is_rule_with_lhs(const Expr* e, const char* lhs_symbol);
bool head_threads_over(const char* h);

/* ------------------------------------------------------------------ */
/* simp_assume.c -- cross-module helpers used by several rewriters     */
/* ------------------------------------------------------------------ */

int  numeric_sign(const Expr* e);
bool fact_is_function(const Expr* f, const char* head, size_t arity);
bool fact_in_domain(const Expr* f, const Expr* x, const char* dom);
bool is_positive_constant_symbol(const char* s);
bool is_real_constant_symbol(const char* s);

bool prov_pos (const AssumeCtx* ctx, const Expr* x);
bool prov_nn  (const AssumeCtx* ctx, const Expr* x);
bool prov_neg (const AssumeCtx* ctx, const Expr* x);
bool prov_np  (const AssumeCtx* ctx, const Expr* x);
bool prov_int (const AssumeCtx* ctx, const Expr* x);
bool prov_re  (const AssumeCtx* ctx, const Expr* x);
bool prov_even(const AssumeCtx* ctx, const Expr* x);

/* ------------------------------------------------------------------ */
/* simp_assume_rewrite.c                                              */
/* ------------------------------------------------------------------ */

Expr* apply_assumption_rules(const Expr* input, const AssumeCtx* ctx);

/* ------------------------------------------------------------------ */
/* simp_trig_roundtrip.c                                              */
/* ------------------------------------------------------------------ */

Expr* transform_trig_roundtrip(const Expr* e);
void  simp_install_roots_of_unity_helpers(void);
Expr* simp_roots_of_unity(const Expr* e);

/* PythagSquareComplete and HalfAngle are defined alongside the
 * trig/exp roundtrip code (they share the memoised seed pattern). */
Expr* transform_pythag_square_complete(const Expr* e);
Expr* transform_halfangle(const Expr* e);

/* ------------------------------------------------------------------ */
/* simp_denest.c                                                      */
/* ------------------------------------------------------------------ */

Expr* simp_radicals(const Expr* e);
Expr* simp_denest_sqrt(const Expr* e, const AssumeCtx* ctx);

/* Defined in simp_denest.c; used cross-module by simp_cuberoot.c. */
bool is_sqrt(const Expr* e);

/* ------------------------------------------------------------------ */
/* simp_cuberoot.c                                                    */
/* ------------------------------------------------------------------ */

Expr* simp_cuberoot(const Expr* e, const AssumeCtx* ctx);

/* ------------------------------------------------------------------ */
/* simp_rationalize.c                                                 */
/* ------------------------------------------------------------------ */

Expr* simp_rationalize_denom(const Expr* e, const AssumeCtx* ctx);

/* ------------------------------------------------------------------ */
/* simp_algebraic.c                                                   */
/* ------------------------------------------------------------------ */

Expr* simp_algebraic(const Expr* e);

/* Used cross-module by simp_builtins.c. */
bool contains_explicit_complex(const Expr* e);

/* ------------------------------------------------------------------ */
/* simp_canon.c                                                       */
/* ------------------------------------------------------------------ */

Expr* simp_lift_common_factor(const Expr* e);
Expr* canon_negate_pairs(const Expr* e);

/* Memoised wrapper used by transforms in many files. Defined in
 * simp_canon.c (the central definition; the forward declaration in
 * simp_denest.c was a duplicate). */
Expr* simp_memo_wrap(const Expr* e, const char* pseudo_head,
                     Expr* (*impl)(const Expr*));

/* ------------------------------------------------------------------ */
/* simp_trig_pi.c                                                     */
/* ------------------------------------------------------------------ */

Expr* simp_trig_pi_canon(const Expr* e);
Expr* transform_pythag_reduce(const Expr* e);
Expr* transform_pythag_canon(const Expr* e);

/* Used by simp_search.c (and the trig_roundtrip Pythag transforms in
 * simp_trig_roundtrip.c) to gate Pythag rewrites. Defined alongside
 * the rest of the Pythag transforms in simp_trig_pi.c. */
bool has_pythag_head(const Expr* e);

/* ------------------------------------------------------------------ */
/* simp_power.c                                                       */
/* ------------------------------------------------------------------ */

Expr* transform_prime_rebase(const Expr* e);
Expr* transform_power_oneify(const Expr* e);
Expr* transform_power_distribute(const Expr* e, const AssumeCtx* ctx);
Expr* transform_radical_canon(const Expr* e);

/* ------------------------------------------------------------------ */
/* simp_tan_add.c                                                     */
/* ------------------------------------------------------------------ */

Expr* transform_tan_addition(const Expr* e);

/* ------------------------------------------------------------------ */
/* simp_logexp_abs.c                                                  */
/* ------------------------------------------------------------------ */

Expr* apply_logexp_rules(const Expr* input, const AssumeCtx* ctx);
Expr* apply_abs_rules(const Expr* input, const AssumeCtx* ctx);
Expr* apply_sqrt_of_square_rules(const Expr* input, const AssumeCtx* ctx);

/* Tree-shape probes used by the simp_search shape classifier. */
bool   contains_abs(const Expr* e);
bool   contains_log(const Expr* e);
bool   contains_power(const Expr* e);
bool   contains_plus_or_times(const Expr* e);
bool   contains_variable(const Expr* e);
bool   contains_trig_or_hyperbolic(const Expr* e);
bool   contains_exp_form(const Expr* e);
bool   contains_sqrt_of_square(const Expr* e);
size_t expr_variables_count_capped(const Expr* e, size_t cap);
bool   ctx_has_facts(const AssumeCtx* ctx);

/* ------------------------------------------------------------------ */
/* simp_search.c                                                      */
/* ------------------------------------------------------------------ */

/* Coarse shape classifier output. Used by builtin_simplify (which takes
 * a SIMP_SHAPE_RATIONAL fast-path through simp_dispatch) and by the
 * specialised pipelines inside simp_search.c. */
typedef enum {
    SIMP_SHAPE_POLYNOMIAL,
    SIMP_SHAPE_RATIONAL,
    SIMP_SHAPE_TRIG,
    SIMP_SHAPE_LOGEXP,
    SIMP_SHAPE_GENERAL
} SimpShape;

SimpShape simp_classify(const Expr* e);

Expr* simp_search(const Expr* original_input, const AssumeCtx* ctx,
                  const Expr* complexity_func);
Expr* simp_dispatch(const Expr* input, const AssumeCtx* ctx,
                    const Expr* complexity_func);

/* Used by simp_search.c gating logic; also called from simp_bottomup.c. */
bool has_non_integer_power(const Expr* e);

/* transform_can_fire is the per-input gate used by both the heuristic
 * search and the bottom-up driver. */
bool transform_can_fire(const char* name, const Expr* e,
                        const AssumeCtx* ctx);

/* ------------------------------------------------------------------ */
/* simp_factorial.c                                                   */
/* ------------------------------------------------------------------ */

Expr* simp_factorial(const Expr* e);
bool  contains_factorial(const Expr* e);

/* ------------------------------------------------------------------ */
/* simp_bottomup.c                                                    */
/* ------------------------------------------------------------------ */

void        simp_memo_init(SimpMemo* m);
void        simp_memo_free(SimpMemo* m);
const Expr* simp_memo_get(SimpMemo* m, const Expr* key);
void        simp_memo_put(SimpMemo* m, const Expr* key, const Expr* value);

Expr* simp_bottomup(const Expr* input, const AssumeCtx* ctx,
                    const Expr* complexity_func, SimpMemo* memo,
                    int depth);

/* ------------------------------------------------------------------ */
/* simp_builtins.c                                                    */
/* ------------------------------------------------------------------ */

/* simp_eq_head_sym is defined here but used cross-module by
 * simp_search.c. */
bool simp_eq_head_sym(const Expr* e, const char* name);

/* is_rational_literal lives in simp_builtins.c and is consumed by many
 * of the rewriters (simp_denest.c, simp_cuberoot.c, simp_rationalize.c,
 * simp_canon.c, simp_search.c). */
bool is_rational_literal(const Expr* e);

/* Public builtins: registered from simp_init. */
Expr* builtin_simplify(Expr* res);
Expr* builtin_simplify_count(Expr* res);
Expr* builtin_assuming(Expr* res);
Expr* builtin_element(Expr* res);

#endif /* MATHILDA_SIMP_INTERNAL_H */
