/* risch_util.h — shared low-level helpers for the transcendental Risch modules.
 *
 * The small evaluation / predicate / polynomial primitives used across the
 * transcendental Risch subsystem (integrate_risch_transcendental.c, its tower
 * and RDE layers, and the case modules split out of it).  Each is a thin,
 * host-grounded wrapper over Mathilda's Expr/eval/poly machinery — build a
 * call, evaluate it, read back a predicate or coefficient.  Defined in
 * risch_util.c.
 *
 * These carry the historical `rt_` prefix ("Risch transcendental") and are the
 * one implementation shared by every module in the subsystem; risch_tower.h and
 * integrate_risch_rde.h include this header so their consumers see the helpers.
 */

#ifndef MATHILDA_RISCH_UTIL_H
#define MATHILDA_RISCH_UTIL_H

#include <stdbool.h>
#include <stddef.h>
#include "expr.h"

/* Build-and-evaluate. */
Expr* rt_eval_call(const char* head, Expr** args, size_t n);  /* adopts args elems */
Expr* rt_eval_own(Expr* e);                                   /* frees e, returns eval */
Expr* rt_eval1(const char* head, Expr* a);
Expr* rt_eval2(const char* head, Expr* a, Expr* b);
Expr* rt_eval3(const char* head, Expr* a, Expr* b, Expr* c);
Expr* rt_template(const char* tmpl, const char** names, Expr** vals, size_t n);

/* Predicates. */
bool  rt_head_is(const Expr* e, const char* name);   /* e == name[...] ? */
bool  rt_is_true(const Expr* e);
bool  rt_is_rat_const(const Expr* e);                /* nonzero rational-number constant */
bool  rt_free_of_x(Expr* e, Expr* x);                /* FreeQ[e, x] */
bool  rt_free_of_head(Expr* e, const char* head);    /* FreeQ[e, head] */
bool  rt_is_poly(Expr* e, Expr* x);                  /* PolynomialQ[e, x] */
bool  rt_is_zero(Expr* e);                           /* exact structural zero test */
bool  rt_has_algebraic_of_x(Expr* e, Expr* x);       /* radical / Surd / Root of x anywhere */

/* Polynomial accessors. */
long  rt_degree(Expr* e, Expr* x);                   /* deg_x, or -1 */
Expr* rt_coeff(Expr* e, Expr* x, long k);            /* Coefficient[e, x, k] */

/* Commensurability-class primitive exponent (see risch_util.c). */
Expr* rt_class_primitive(Expr** ws, size_t nw, long* kof);

/* First x-dependent kernel argument/exponent (borrowed), or NULL. */
Expr* rt_find_log_of_x(Expr* e, Expr* x);
Expr* rt_find_exp_of_x(Expr* e, Expr* x);

/* True iff kernel-defining function u is a rational function of x alone
 * (no nested exp/log of x) — i.e. a valid single extension. */
bool  rt_kernel_simple(Expr* u, Expr* x);

/* Diff-back verifier: is Simplify[D[result,x] - f] identically 0? */
bool  rt_verify_antideriv(Expr* result, Expr* f, Expr* x);

#endif /* MATHILDA_RISCH_UTIL_H */
