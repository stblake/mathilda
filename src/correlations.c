/* Mathilda — ListCorrelate.
 *
 * A thin front end: ListCorrelate shares the general engine in convolutions.c,
 * differing from ListConvolve only in the alignment direction (Sum_r K_r a_{s+r}
 * vs Sum_r K_r a_{s-r}) and in the default/negated overhang conventions, all of
 * which the engine handles via its ConvMode parameter. See convolutions.h. */

#include "correlations.h"
#include "convolutions.h"
#include "symtab.h"
#include "attr.h"

Expr* builtin_list_correlate(Expr* res) {
    return conv_engine(res, CONV_MODE_CORRELATE);
}

void correlations_init(void) {
    symtab_add_builtin("ListCorrelate", builtin_list_correlate);
    symtab_get_def("ListCorrelate")->attributes |= ATTR_PROTECTED;
}
