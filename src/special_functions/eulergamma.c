/* Mathilda -- EulerGamma, the Euler-Mascheroni constant gamma.
 *
 *   EulerGamma = lim_{n->oo} (HarmonicNumber[n] - Log[n])
 *              ~= 0.5772156649015328606065120900824024310421593359399...
 *
 * Like Pi and E, EulerGamma is propagated as an exact, unevaluated symbol;
 * its numeric value is materialised only on demand by N[] (machine or
 * arbitrary precision) -- see the constant table in src/numeric.c, which
 * routes the MPFR path through mpfr_const_euler for full-precision digits.
 *
 * This module owns nothing more than the symbol's *identity*: it marks the
 * symbol Constant (so D treats it as a constant and it reads as a genuine
 * mathematical constant) and Protected (so it cannot be reassigned). All of
 * the actual numeric / derivative / NumericQ behaviour is supplied by the
 * generic subsystems documented in eulergamma.h. */

#include "eulergamma.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"

void eulergamma_init(void) {
    /* Ensure a SymbolDef exists, then stamp the Mathematica attributes.
     * EulerGamma has no DownValues / builtin function pointer -- it is a
     * pure constant symbol. */
    SymbolDef* def = symtab_get_def(SYM_EulerGamma);
    if (def) {
        def->attributes |= (ATTR_CONSTANT | ATTR_PROTECTED);
    }
    /* Docstring lives in src/info.c (info_init). */
}
