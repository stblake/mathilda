/* Mathilda -- StieltjesGamma, the Stieltjes constants gamma_n.
 *
 *   StieltjesGamma[n] = gamma_n, the coefficients of the Laurent expansion
 *   of the Riemann zeta function about s = 1:
 *
 *     zeta(s) = 1/(s-1) + Sum_{n>=0} ((-1)^n / n!) gamma_n (s-1)^n.
 *
 * gamma_0 is the Euler-Mascheroni constant EulerGamma. The higher constants
 * have no elementary closed form, so StieltjesGamma is inert: it stays
 * symbolic, except for the single reduction StieltjesGamma[0] -> EulerGamma.
 * It is the natural output of Series[Zeta[x], {x, 1, n}] (and the Taylor
 * expansion of Zeta about 0).
 *
 * Like LogGamma (see src/polygamma.c), this module owns only the symbol's
 * identity and the n = 0 reduction; all generic symbolic behaviour comes from
 * the evaluator.
 *
 * Attributes: Listable, Protected. */
#include "stieltjesgamma.h"

#include <stdbool.h>

#include "attr.h"
#include "symtab.h"

/* StieltjesGamma[0] -> EulerGamma; everything else stays symbolic. */
Expr* builtin_stieltjesgamma(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* n = res->data.function.args[0];
    if (n->type == EXPR_INTEGER && n->data.integer == 0)
        return expr_new_symbol("EulerGamma");
    return NULL; /* inert */
}

void stieltjesgamma_init(void) {
    symtab_add_builtin("StieltjesGamma", builtin_stieltjesgamma);
    symtab_get_def("StieltjesGamma")->attributes |=
        (ATTR_LISTABLE | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
