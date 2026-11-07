/* RotationMatrix -- 2-D and 3-D rotation matrices.
 *
 *   RotationMatrix[t]              {{Cos t, -Sin t}, {Sin t, Cos t}}
 *   RotationMatrix[{u, v}]         2-D or 3-D: the rotation taking the
 *                                  direction of u to that of v, in the plane
 *                                  they span (by the angle between them)
 *   RotationMatrix[t, w]           3-D: rotation by t about the axis w
 *                                  (counter-clockwise looking down w)
 *   RotationMatrix[t, {u, v}]      rotation by t in the plane of u and v, in
 *                                  the sense that turns u toward v
 *
 * Every 3-D form is Rodrigues' formula with the axis left UNNORMALISED, so an
 * exact input gives an exact matrix without a radical it does not need:
 *
 *   R = c I + (s / |w|) [w]x + ((1 - c) / w.w) w w^T,
 *
 * and for {u, v}: w = u x v, c = u.v / N, s / |w| = 1 / N, N = Sqrt[(u.u)(v.v)]
 * -- RotationMatrix[{{1,2,3},{3,2,1}}] is all rationals (16/21, ...).
 *
 * The entries are built as expressions and evaluated, so symbolic angles and
 * vectors work.  Deliberate difference: a symbolic axis is taken to be real
 * (the reference's general form carries Conjugate[] of each component); for
 * real input the two agree.  Only 2-D and 3-D are supported (n-D stays
 * unevaluated).
 *
 * Messages, as in the reference: RotationMatrix::inpf (not a scalar or a pair
 * of vectors), ::v3d (an axis not of length 3), ::spln (u, v do not define a
 * plane: anti-parallel 3-D vectors).  No numeric fast path: the result is a
 * small matrix, not an element-wise numeric map (documented EXEMPT).
 *
 * Attributes: Protected.
 */

#include "linalg.h"
#include "common.h"       /* builtin_arg_error */
#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "message.h"
#include "print.h"
#include "sym_names.h"
#include "symtab.h"
#include <stdbool.h>
#include <stdlib.h>

/* ---- expression helpers (each consumes its arguments) ----------------- */

static Expr* rm_fn(const char* head, Expr** a, size_t n) {
    return expr_new_function(expr_new_symbol(head), a, n);
}
static Expr* rm_add(Expr* a, Expr* b) { return rm_fn(SYM_Plus, (Expr*[]){ a, b }, 2); }
static Expr* rm_mul(Expr* a, Expr* b) { return rm_fn(SYM_Times, (Expr*[]){ a, b }, 2); }
static Expr* rm_neg(Expr* a) { return rm_mul(expr_new_integer(-1), a); }
static Expr* rm_pow(Expr* a, Expr* b) { return rm_fn(SYM_Power, (Expr*[]){ a, b }, 2); }
static Expr* rm_div(Expr* a, Expr* b) { return rm_mul(a, rm_pow(b, expr_new_integer(-1))); }
static Expr* rm_sqrt(Expr* a) {
    return rm_pow(a, rm_fn(SYM_Rational, (Expr*[]){ expr_new_integer(1), expr_new_integer(2) }, 2));
}
static Expr* rm_c(const Expr* e) { return expr_copy((Expr*)e); }

static bool rm_is_list(const Expr* e, size_t n) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List &&
           (n == 0 || e->data.function.arg_count == n);
}
static bool rm_is_vector(const Expr* e) {
    if (!rm_is_list(e, 0) || e->data.function.arg_count == 0) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rm_is_list(e->data.function.args[i], 0)) return false;
    return true;
}
static const Expr* rm_at(const Expr* v, size_t i) { return v->data.function.args[i]; }

/* u . v, unevaluated. */
static Expr* rm_dot(const Expr* u, const Expr* v) {
    size_t n = u->data.function.arg_count;
    Expr** t = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) t[i] = rm_mul(rm_c(rm_at(u, i)), rm_c(rm_at(v, i)));
    Expr* r = rm_fn(SYM_Plus, t, n);
    free(t);
    return r;
}

/* True iff the evaluated `e` is an exact or inexact numeric zero. */
static bool rm_is_zero(const Expr* e) {
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_REAL) return e->data.real == 0.0;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return mpfr_zero_p(e->data.mpfr) != 0;
#endif
    return false;
}

/* True iff the evaluated `e` is a negative real number. */
static bool rm_is_negative(const Expr* e) {
    Expr* t = eval_and_free(rm_fn(SYM_Negative, (Expr*[]){ rm_c(e) }, 1));
    bool r = t->type == EXPR_SYMBOL && t->data.symbol.name == SYM_True;
    expr_free(t);
    return r;
}

/* {{c, -s}, {s, c}}; consumes c and s. */
static Expr* rm_2d(Expr* c, Expr* s) {
    Expr* r0 = rm_fn(SYM_List, (Expr*[]){ rm_c(c), rm_neg(rm_c(s)) }, 2);
    Expr* r1 = rm_fn(SYM_List, (Expr*[]){ s, c }, 2);
    return eval_and_free(rm_fn(SYM_List, (Expr*[]){ r0, r1 }, 2));
}

/* Rodrigues: c I + k [w]x + m w w^T, with k = s/|w| and m = (1-c)/w.w given;
 * consumes c, k, m. */
static Expr* rm_3d(const Expr* w, Expr* c, Expr* k, Expr* m) {
    const Expr* w1 = rm_at(w, 0); const Expr* w2 = rm_at(w, 1); const Expr* w3 = rm_at(w, 2);
    /* [w]x = {{0, -w3, w2}, {w3, 0, -w1}, {-w2, w1, 0}} */
    Expr* cross[3][3] = {
        { expr_new_integer(0), rm_neg(rm_c(w3)), rm_c(w2) },
        { rm_c(w3), expr_new_integer(0), rm_neg(rm_c(w1)) },
        { rm_neg(rm_c(w2)), rm_c(w1), expr_new_integer(0) },
    };
    Expr* rows[3];
    for (int i = 0; i < 3; i++) {
        Expr* row[3];
        for (int j = 0; j < 3; j++) {
            Expr* e = rm_add(rm_mul(rm_c(k), cross[i][j]),
                             rm_mul(rm_c(m), rm_mul(rm_c(rm_at(w, (size_t)i)), rm_c(rm_at(w, (size_t)j)))));
            if (i == j) e = rm_add(rm_c(c), e);
            row[j] = e;
        }
        rows[i] = rm_fn(SYM_List, row, 3);
    }
    expr_free(c); expr_free(k); expr_free(m);
    return eval_and_free(rm_fn(SYM_List, rows, 3));
}

/* w = u x v (3-D), unevaluated components in a List. */
static Expr* rm_cross(const Expr* u, const Expr* v) {
    Expr* comp[3];
    for (int i = 0; i < 3; i++) {
        size_t a = (size_t)((i + 1) % 3), b = (size_t)((i + 2) % 3);
        comp[i] = rm_add(rm_mul(rm_c(rm_at(u, a)), rm_c(rm_at(v, b))),
                         rm_neg(rm_mul(rm_c(rm_at(u, b)), rm_c(rm_at(v, a)))));
    }
    return eval_and_free(rm_fn(SYM_List, comp, 3));
}

static void rm_msg(const char* tag, const char* fmt, const Expr* a, const Expr* b) {
    char* sa = a ? expr_to_string((Expr*)a) : NULL;
    char* sb = b ? expr_to_string((Expr*)b) : NULL;
    mth_message("RotationMatrix", tag, fmt, sa ? sa : "", sb ? sb : "");
    free(sa); free(sb);
}

/* The rotation taking direction u to direction v (2-D or 3-D). */
static Expr* rm_two_vectors(const Expr* u, const Expr* v) {
    size_t n = u->data.function.arg_count;
    Expr* uu = eval_and_free(rm_dot(u, u));
    Expr* vv = eval_and_free(rm_dot(v, v));
    if (rm_is_zero(uu) || rm_is_zero(vv)) { expr_free(uu); expr_free(vv); return NULL; }
    Expr* N = eval_and_free(rm_sqrt(rm_mul(uu, vv)));      /* |u| |v| */
    Expr* c = eval_and_free(rm_div(rm_dot(u, v), rm_c(N)));
    if (n == 2) {
        Expr* x = rm_add(rm_mul(rm_c(rm_at(u, 0)), rm_c(rm_at(v, 1))),
                         rm_neg(rm_mul(rm_c(rm_at(u, 1)), rm_c(rm_at(v, 0)))));
        Expr* s = eval_and_free(rm_div(x, N));
        return rm_2d(c, s);
    }
    Expr* w = rm_cross(u, v);
    Expr* ww = eval_and_free(rm_dot(w, w));
    if (rm_is_zero(ww)) {
        /* parallel: the identity; anti-parallel: no plane is defined */
        bool anti = rm_is_negative(c);
        expr_free(w); expr_free(ww); expr_free(N); expr_free(c);
        if (anti) {
            rm_msg("spln", "Vectors %s and %s do not define a plane.", u, v);
            return NULL;
        }
        return eval_and_free(rm_fn("IdentityMatrix", (Expr*[]){ expr_new_integer(3) }, 1));
    }
    Expr* k = eval_and_free(rm_pow(N, expr_new_integer(-1)));
    Expr* m = eval_and_free(rm_div(rm_add(expr_new_integer(1), rm_neg(rm_c(c))), ww));
    Expr* r = rm_3d(w, c, k, m);
    expr_free(w);
    return r;
}

/* Rotation by t about the axis w (3-D). */
static Expr* rm_axis(const Expr* t, const Expr* w) {
    Expr* ww = eval_and_free(rm_dot(w, w));
    if (rm_is_zero(ww)) { expr_free(ww); return NULL; }
    Expr* c = eval_and_free(rm_fn(SYM_Cos, (Expr*[]){ rm_c(t) }, 1));
    Expr* s = rm_fn(SYM_Sin, (Expr*[]){ rm_c(t) }, 1);
    Expr* k = eval_and_free(rm_div(s, rm_sqrt(rm_c(ww))));
    Expr* m = eval_and_free(rm_div(rm_add(expr_new_integer(1), rm_neg(rm_c(c))), ww));
    return rm_3d(w, c, k, m);
}

Expr* builtin_rotationmatrix(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** a = res->data.function.args;

    if (argc == 1) {
        Expr* x = a[0];
        if (rm_is_list(x, 2) && rm_is_vector(rm_at(x, 0)) && rm_is_vector(rm_at(x, 1))) {
            const Expr* u = rm_at(x, 0); const Expr* v = rm_at(x, 1);
            size_t n = u->data.function.arg_count;
            if (v->data.function.arg_count != n || (n != 2 && n != 3)) return NULL;
            return rm_two_vectors(u, v);
        }
        if (rm_is_list(x, 0)) {
            rm_msg("inpf", "%s is not a scalar or a list of two vectors.", x, NULL);
            return NULL;
        }
        Expr* c = eval_and_free(rm_fn(SYM_Cos, (Expr*[]){ rm_c(x) }, 1));
        Expr* s = eval_and_free(rm_fn(SYM_Sin, (Expr*[]){ rm_c(x) }, 1));
        return rm_2d(c, s);
    }
    if (argc == 2) {
        Expr* t = a[0];
        Expr* w = a[1];
        if (rm_is_list(t, 0)) return NULL;
        if (rm_is_list(w, 2) && rm_is_vector(rm_at(w, 0)) && rm_is_vector(rm_at(w, 1))) {
            /* the plane of u and v, turning u toward v */
            const Expr* u = rm_at(w, 0); const Expr* v = rm_at(w, 1);
            size_t n = u->data.function.arg_count;
            if (v->data.function.arg_count != n || (n != 2 && n != 3)) return NULL;
            if (n == 2) {
                Expr* x = eval_and_free(rm_add(rm_mul(rm_c(rm_at(u, 0)), rm_c(rm_at(v, 1))),
                                     rm_neg(rm_mul(rm_c(rm_at(u, 1)), rm_c(rm_at(v, 0))))));
                bool neg = rm_is_negative(x);
                bool zero = rm_is_zero(x);
                expr_free(x);
                if (zero) { rm_msg("spln", "Vectors %s and %s do not define a plane.", u, v); return NULL; }
                Expr* tt = neg ? rm_neg(rm_c(t)) : rm_c(t);
                Expr* c = eval_and_free(rm_fn(SYM_Cos, (Expr*[]){ rm_c(tt) }, 1));
                Expr* s = eval_and_free(rm_fn(SYM_Sin, (Expr*[]){ tt }, 1));
                return rm_2d(c, s);
            }
            Expr* ax = rm_cross(u, v);
            Expr* ww = eval_and_free(rm_dot(ax, ax));
            bool zero = rm_is_zero(ww);
            expr_free(ww);
            if (zero) {
                rm_msg("spln", "Vectors %s and %s do not define a plane.", u, v);
                expr_free(ax);
                return NULL;
            }
            Expr* r = rm_axis(t, ax);
            expr_free(ax);
            return r;
        }
        if (rm_is_vector(w)) {
            if (w->data.function.arg_count != 3) {
                rm_msg("v3d", "%s is not a vector of length 3.", w, NULL);
                return NULL;
            }
            return rm_axis(t, w);
        }
        return NULL;
    }
    return builtin_arg_error("RotationMatrix", argc, 1, 2);
}

void rotationmatrix_init(void) {
    symtab_add_builtin("RotationMatrix", builtin_rotationmatrix);
    symtab_get_def("RotationMatrix")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("RotationMatrix",
        "RotationMatrix[t]\n"
        "\tthe 2-D matrix {{Cos[t], -Sin[t]}, {Sin[t], Cos[t]}} rotating by t radians.\n"
        "RotationMatrix[t, w]\n"
        "\tthe 3-D rotation by t counter-clockwise about the axis w.\n"
        "RotationMatrix[{u, v}]\n"
        "\tthe 2-D or 3-D rotation taking the direction of u to that of v.\n"
        "RotationMatrix[t, {u, v}]\n"
        "\tthe rotation by t in the plane of u and v, turning u toward v.\n"
        "Exact input gives exact entries (Rodrigues' formula with the axis unnormalised).");
}
