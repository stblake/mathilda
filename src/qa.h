/* qa.h — Q(α) algebraic-number type (Phase G1, Trager's algorithm).
 *
 * QAExt is an algebraic extension Q(α) = Q[α] / (P_α(α)) defined by a
 * minimal polynomial P_α ∈ Q[α]. The minimal polynomial must be
 * irreducible over Q for inversion to always succeed; we don't enforce
 * irreducibility here (the caller is responsible for that contract).
 *
 * QANum represents an element of Q(α), stored as a polynomial in α of
 * degree < deg P_α. Arithmetic operations preserve reduced form by
 * reducing modulo P_α whenever a multiplication would otherwise raise
 * the degree.
 *
 * Storage convention: coef is an mpq_t array of exactly ext->deg
 * entries; coef[i] is the coefficient of α^i. Lower-degree elements
 * are zero-padded out to ext->deg. The QAExt pointer is borrowed —
 * the caller owns the extension throughout the QANum's lifetime.
 *
 * Phase G1 is the substrate for Trager's algebraic-factoring algorithm
 * (Phase G in FACTOR_PLAN.md §14). Phase G2 builds Q(α)[x] on top of
 * it; Phase G3 computes norms via Resultant_y(P_α, f). */

#ifndef MATHILDA_QA_H
#define MATHILDA_QA_H

#include <stddef.h>
#include <stdbool.h>
#include <gmp.h>

struct Expr;

/* The extension Q(α) defined by a minimal polynomial of degree `deg`.
 * coef[i] is the coefficient of y^i in P_α(y), 0 ≤ i ≤ deg. */
typedef struct QAExt {
    size_t deg;
    mpq_t* coef;
} QAExt;

/* An element of Q(α). coef[i] is the coefficient of α^i, 0 ≤ i < ext->deg. */
typedef struct QANum {
    const QAExt* ext;
    mpq_t* coef;
} QANum;

/* ===================== QAExt construction / destruction ================== */

/* Allocate a QAExt of the given degree. All coefficients start at 0;
 * caller must populate them via qaext_set_coef_si / qaext_set_coef_mpq
 * before using any QANum operation that reduces modulo P_α. */
QAExt* qaext_new(size_t deg);

/* Release all storage owned by `ext`, including the mpq_t coefficients. */
void qaext_free(QAExt* ext);

/* Populate the i-th coefficient of P_α from a rational num/den or mpq_t. */
void qaext_set_coef_si(QAExt* ext, size_t i, long num, long den);
void qaext_set_coef_mpq(QAExt* ext, size_t i, const mpq_t v);

/* Convenience constructors for the two extensions Phase G5 cares about:
 *
 *   qaext_sqrt_si(c)         P_α(y) = y^2 - c   (degree 2)
 *   qaext_root_si(c, n)      P_α(y) = y^n - c   (degree n)
 *
 * Returns a fresh QAExt with the leading coefficient set to 1. */
QAExt* qaext_sqrt_si(long c);
QAExt* qaext_root_si(long c, unsigned n);

/* =========================== QANum construction ========================== */

/* Construct elementary values. `ext` is borrowed (must outlive the QANum). */
QANum* qa_zero(const QAExt* ext);
QANum* qa_one(const QAExt* ext);
QANum* qa_alpha(const QAExt* ext);                 /* α itself */
QANum* qa_from_mpq(const QAExt* ext, const mpq_t v);
QANum* qa_from_si(const QAExt* ext, long num, long den);

/* α^p ∈ Q(α) for any integer p, reduced modulo P_α.  For p < 0 uses
 * `qa_inverse` (returns NULL iff `qa_inverse` fails — i.e. P_α is
 * reducible and α^|p| is a zero divisor).  Repeated-squaring; O(log |p|)
 * multiplications. */
QANum* qa_alpha_power_signed(const QAExt* ext, long p);

/* Deep copy. Shares the extension by reference. */
QANum* qa_copy(const QANum* a);

void   qa_free(QANum* a);

/* ============================== Predicates ============================== */

bool qa_is_zero(const QANum* a);
bool qa_is_one(const QANum* a);
bool qa_eq(const QANum* a, const QANum* b);        /* assumes same ext */

/* ============================ Arithmetic ============================== */
/* All operations return a freshly-allocated QANum and leave their
 * arguments unmodified. They share the same extension; mixing
 * QANums from different extensions is undefined. */

QANum* qa_add(const QANum* a, const QANum* b);
QANum* qa_sub(const QANum* a, const QANum* b);
QANum* qa_neg(const QANum* a);
QANum* qa_mul(const QANum* a, const QANum* b);     /* reduces mod P_α */

QANum* qa_scale_mpq(const QANum* a, const mpq_t k);
QANum* qa_scale_si(const QANum* a, long num, long den);

/* Multiplicative inverse via extended Euclidean over Q[y].
 * Returns NULL if a == 0 or if gcd(a, P_α) != 1 (which can only
 * happen when P_α is reducible and a is a zero divisor). */
QANum* qa_inverse(const QANum* a);

/* a / b  =  a * inverse(b). Returns NULL if b is zero / non-invertible. */
QANum* qa_div(const QANum* a, const QANum* b);

#endif
