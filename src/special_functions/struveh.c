/* Mathilda -- the Struve function H_nu(z).
 *
 *   StruveH[nu, z]   H_nu(z) = (z/2)^(nu+1) (2/(Sqrt[Pi] Gamma[nu+3/2]))
 *                              1F2(1; 3/2, nu+3/2; -z^2/4)
 *
 * Evaluation:
 *   StruveH[nu, 0] = 0            exact, for a real nu > -1
 *   numeric nu and inexact z      via the 1F2 representation (reusing the tested
 *                                 HypergeometricPFQ numeric evaluator)
 *   everything else               stays symbolic (return NULL)
 *
 * The 1F2 identity is also the reduction used by the half-line Mellin transform
 * (reduce_to_hypergeometric in integrate_ramanujan.c), so Integrate[x^(s-1)
 * StruveH[0, a x], {x,0,Infinity}] closes through rec_pfq.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "struveh.h"

#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "parse.h"
#include "sym_names.h"
#include "common.h"
#include "numeric.h"

#include <stdbool.h>

/* Small expr builders (local idiom). */
static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_fn1(const char* name, Expr* a) {
    return expr_new_function(expr_new_symbol(name), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* name, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(name), (Expr*[]){ a, b }, 2);
}
/* evaluate `call` to a fixed point, then free it (evaluate borrows its arg). */
static Expr* sh_eval(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* H_nu(z) via its 1F2 representation, numericised; NULL if the result is not a
 * number (e.g. the pFq evaluator declined). */
static Expr* struveh_numeric(const Expr* nu, const Expr* z) {
    static const char* TMPL =
        "(zz/2)^(nn+1) (2/(Sqrt[Pi] Gamma[nn+3/2])) "
        "HypergeometricPFQ[{1}, {3/2, nn+3/2}, -zz^2/4]";
    Expr* tmpl = parse_expression(TMPL);
    if (!tmpl) return NULL;
    Expr* rules = expr_new_function(expr_new_symbol("List"), (Expr*[]){
        mk_fn2("Rule", expr_new_symbol("nn"), expr_copy((Expr*)nu)),
        mk_fn2("Rule", expr_new_symbol("zz"), expr_copy((Expr*)z)) }, 2);
    Expr* subbed = sh_eval(mk_fn2("ReplaceAll", tmpl, rules));
    Expr* r = sh_eval(mk_fn1("N", subbed));
    if (r && expr_is_numeric_quantity(r)) return r;
    if (r) expr_free(r);
    return NULL;
}

Expr* builtin_struveh(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* nu = res->data.function.args[0];
    Expr* z  = res->data.function.args[1];

    /* H_nu(0) = 0 for a real nu > -1 (exact). */
    if (z->type == EXPR_INTEGER && z->data.integer == 0) {
        double nuv;
        if (common_machine_real_value(nu, &nuv) && nuv > -1.0) return mk_int(0);
    }

    /* Numeric: a numeric order nu and an inexact (N-triggered) real z. */
    if (z->type == EXPR_REAL) {
        double nuv;
        if (common_machine_real_value(nu, &nuv)) {
            Expr* v = struveh_numeric(nu, z);
            if (v) return v;
        }
    }
    return NULL;
}

void struveh_init(void) {
    symtab_add_builtin("StruveH", builtin_struveh);
    symtab_get_def("StruveH")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
