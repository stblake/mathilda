/* Mathilda — Compile[] engine: the LinkedList opcodes at run time.
 *
 * Each op works on the list's typed storage through the same ll_* helpers the
 * interpreter uses (src/linkedlist.h), so a compiled push, read or Total is the
 * interpreted one by construction: the same slots, the same kind rules, the same
 * summation tree. No Expr is built per element.
 *
 * A list register holds an EXPR_LINKEDLIST node in `.arr`. The emitter
 * guarantees its storage kind matches the op's element type (LLF_ELEM): a list
 * the body created only ever receives pushes of that type, and an argument or a
 * folded global was checked at the boundary. A list that still turns out to be
 * of another kind fails the op, and the interpreter answers.
 *
 * Every failure returns false: the call is abandoned and the caller's undo frame
 * reverts any change to a list the run did not create (compile_vm.c,
 * compiled_eval). That is how Pop[] of an empty list, an index out of range and
 * an Integer Total past int64 reach the interpreter's own answer and message
 * without the compiled run leaving a trace.
 */
#include "compile_internal.h"
#include "../linkedlist.h"

#include <complex.h>
#include <string.h>

static LinkedList* reg_list(const Slot* s) { return expr_linkedlist(s->arr); }

/* The storage kind an element type is held as. */
static LLKind elem_kind(CompileType e) {
    return e == CT_INT ? LL_I64 : e == CT_REAL ? LL_REAL : LL_CPLX;
}

/* A 1-based Part index (negative from the end) to a 0-based slot, or -1. */
static int64_t resolve_index(const LinkedList* ll, long long k) {
    int64_t n = ll_length(ll);
    int64_t j = k > 0 ? (int64_t)k - 1 : n + (int64_t)k;
    return (k == 0 || j < 0 || j >= n) ? -1 : j;
}

static bool read_elem(const LinkedList* ll, int64_t j, CompileType e, Slot* out) {
    switch (e) {
        case CT_INT: { int64_t v; if (!ll_get_i64(ll, j, &v)) return false; out->i = (long long)v; return true; }
        case CT_REAL: return ll_get_real(ll, j, &out->r);
        case CT_COMPLEX: {
            double re, im;
            if (!ll_get_cplx(ll, j, &re, &im)) return false;
            out->z = re + im * I;
            return true;
        }
        default: return false;
    }
}

static void push_elem(LinkedList* ll, bool front, CompileType e, const Slot* v) {
    switch (e) {
        case CT_INT:
            if (front) ll_push_front_i64(ll, (int64_t)v->i); else ll_push_back_i64(ll, (int64_t)v->i);
            break;
        case CT_REAL:
            if (front) ll_push_front_real(ll, v->r); else ll_push_back_real(ll, v->r);
            break;
        default:
            if (front) ll_push_front_cplx(ll, creal(v->z), cimag(v->z));
            else       ll_push_back_cplx(ll, creal(v->z), cimag(v->z));
            break;
    }
}

/* Replace the list in a destination register (a reused register may still hold
 * last call's list). */
static void set_list(Slot* d, LinkedList* ll) {
    expr_free(d->arr);
    d->arr = expr_new_linkedlist(ll);
}

static void free_consumed(unsigned f, Slot* a, Slot* b) {
    if ((f & LLF_FREE_A) && a->arr) { expr_free(a->arr); a->arr = NULL; }
    if ((f & LLF_FREE_B) && b->arr) { expr_free(b->arr); b->arr = NULL; }
}

/* Take[l, k] / Drop[l, k] for one integer k, as a range [lo, hi) to keep; false
 * when k is outside the list (the interpreter then gives Take::take). */
static bool take_drop_range(const LinkedList* src, long long k, bool take,
                            int64_t* lo1, int64_t* hi1, int64_t* lo2, int64_t* hi2) {
    int64_t n = ll_length(src);
    if (k > n || -k > n) return false;
    int64_t a = k >= 0 ? 0 : n + (int64_t)k;        /* positions named: [a, b) */
    int64_t b = k >= 0 ? (int64_t)k : n;
    if (take) { *lo1 = a; *hi1 = b; *lo2 = *hi2 = 0; }
    else      { *lo1 = 0; *hi1 = a; *lo2 = b; *hi2 = n; }
    return true;
}

static bool llist_struct(const Instr* c, Slot* d, Slot* a, Slot* b) {
    const LinkedList* src = reg_list(a);
    if (!src) return false;
    CompileType e = LLF_ELEM(c->flags);
    int64_t n = ll_length(src);
    LinkedList* out = ll_new();
    switch ((int)c->imm.i) {
        case LLS_MOST:
        case LLS_REST:
            if (n == 0) { ll_release(out); return false; }     /* Most::nomost */
            if (c->imm.i == LLS_MOST) ll_append_range(out, src, 0, n - 1);
            else                      ll_append_range(out, src, 1, n);
            break;
        case LLS_REVERSE:
            ll_append_range_rev(out, src, 0, n);
            break;
        case LLS_TAKE:
        case LLS_DROP: {
            int64_t l1, h1, l2, h2;
            if (!take_drop_range(src, b->i, c->imm.i == LLS_TAKE, &l1, &h1, &l2, &h2)) {
                ll_release(out); return false;
            }
            ll_append_range(out, src, l1, h1);
            ll_append_range(out, src, l2, h2);
            break;
        }
        case LLS_APPEND:
        case LLS_PREPEND:
            ll_append_range(out, src, 0, n);
            push_elem(out, c->imm.i == LLS_PREPEND, e, b);
            break;
        case LLS_JOIN: {
            const LinkedList* src2 = reg_list(b);
            if (!src2 || src2 == out) { ll_release(out); return false; }
            ll_append_range(out, src, 0, n);
            ll_append_range(out, src2, 0, ll_length(src2));
            break;
        }
        default:
            ll_release(out);
            return false;
    }
    free_consumed(c->flags, a, b);
    set_list(d, out);
    return true;
}

bool vm_llist_op(const Instr* c, Slot* d, Slot* a, Slot* b) {
    const unsigned f = c->flags;
    const CompileType e = LLF_ELEM(f);
    switch (c->op) {
        case OP_LL_NEW:
            set_list(d, ll_new());
            return true;

        case OP_LL_FROMARR: {
            LinkedList* ll = ll_new();
            if (!a->arr || !ll_append_elements(ll, a->arr)) { ll_release(ll); return false; }
            free_consumed(f, a, b);
            set_list(d, ll);
            return true;
        }

        case OP_LL_RETAIN: {
            if (!reg_list(a)) return false;
            Expr* x = expr_copy(a->arr);
            expr_free(d->arr);
            d->arr = x;
            return true;
        }

        case OP_LL_REF: {
            Expr* x = expr_copy((Expr*)c->imm.p);
            expr_free(d->arr);
            d->arr = x;
            return true;
        }

        case OP_LL_PUSH: {
            LinkedList* ll = reg_list(d);
            if (!ll) return false;
            if (ll_length(ll) > 0 && ll_kind(ll) != elem_kind(e)) return false;
            push_elem(ll, (f & LLF_FRONT) != 0, e, a);
            return true;
        }

        case OP_LL_POP: {
            LinkedList* ll = reg_list(a);
            if (!ll || ll_kind(ll) != elem_kind(e)) return false;   /* empty: Pop::empty */
            double raw[2];
            if (!ll_pop_typed(ll, (f & LLF_FRONT) != 0, raw)) return false;
            if (e == CT_INT) { int64_t v; memcpy(&v, raw, sizeof v); d->i = (long long)v; }
            else if (e == CT_REAL) d->r = raw[0];
            else d->z = raw[0] + raw[1] * I;
            return true;
        }

        case OP_LL_GET: {
            const LinkedList* ll = reg_list(a);
            if (!ll) return false;
            int64_t j = resolve_index(ll, b->i);
            return j >= 0 && read_elem(ll, j, e, d);
        }

        case OP_LL_SET: {
            LinkedList* ll = reg_list(d);
            if (!ll) return false;
            int64_t j = resolve_index(ll, a->i);
            if (j < 0) return false;                                 /* Set::partw */
            switch (e) {
                case CT_INT:  return ll_set_i64(ll, j, (int64_t)b->i);
                case CT_REAL: return ll_set_real(ll, j, b->r);
                default:      return ll_set_cplx(ll, j, creal(b->z), cimag(b->z));
            }
        }

        case OP_LL_LEN: {
            const LinkedList* ll = reg_list(a);
            if (!ll) return false;
            d->i = (long long)ll_length(ll);
            return true;
        }

        case OP_LL_TOTAL: {
            /* The interpreter's own Total: exact int64 (a sum past int64 is a
             * bignum, which no machine register holds: fail), or the packed
             * pairwise tree for reals and complexes, bit for bit. */
            const LinkedList* ll = reg_list(a);
            if (!ll) return false;
            if (ll_length(ll) == 0) {
                if (e == CT_INT) d->i = 0; else if (e == CT_REAL) d->r = 0.0; else d->z = 0.0;
                return e == CT_INT;            /* Total[{}] is the Integer 0 */
            }
            if (ll_kind(ll) != elem_kind(e)) return false;
            Expr* s = ll_sum_machine(ll);
            if (!s) return false;
            bool ok = true;
            if (e == CT_INT) {
                if (s->type == EXPR_INTEGER) d->i = (long long)s->data.integer; else ok = false;
            } else if (e == CT_REAL) {
                if (s->type == EXPR_REAL) d->r = s->data.real; else ok = false;
            } else {
                if (s->type == EXPR_FUNCTION && s->data.function.arg_count == 2
                    && s->data.function.args[0]->type == EXPR_REAL
                    && s->data.function.args[1]->type == EXPR_REAL)
                    d->z = s->data.function.args[0]->data.real
                         + s->data.function.args[1]->data.real * I;
                else ok = false;
            }
            expr_free(s);
            return ok;
        }

        case OP_LL_NORMAL: {
            const LinkedList* ll = reg_list(a);
            if (!ll || ll_kind(ll) != elem_kind(e)) return false;    /* empty: {} */
            Expr* r = ll_to_ndarray(ll);
            if (!r) return false;
            expr_free(d->arr);
            d->arr = r;
            return true;
        }

        case OP_LL_STRUCT:
            return llist_struct(c, d, a, b);

        default:
            return false;
    }
}
