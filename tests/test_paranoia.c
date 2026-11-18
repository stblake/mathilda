/* Mathilda -- W. Kahan's PARANOIA floating-point test, as a ctest gate
 * (numeric-tower review, phase 9; docs/design/numeric_tower_review.md).
 *
 * The test itself is tests/paranoia/paranoia.m, a port of the Sumner/Gay C
 * translation: milestones 0-100 line by line, 110-220 adapted to Mathilda's
 * machine/MPFR number model (see the header of that file). This driver Get[]s
 * it and calls Paranoia`ParanoiaGate[], which runs the modes "Machine",
 * {"MPFR", 30|50|100} and "Mixed" with exact comparisons, prints every report
 * (FAILURE / SERIOUS / DEFECT / FLAW counts and messages, plus the documented
 * policy outcomes), and fails when a FAILURE or SERIOUS appears that is not on
 * the ratcheted known-bug list in paranoia.m. PARANOIA_STRICT=1 ignores that
 * list. The whole run takes about a second.
 *
 * argv[1]: path to paranoia.m (CMake passes the source-tree path).
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* Evaluate `src` and return its printed form (caller frees). */
static char* eval_str(const char* src) {
    Expr* e = parse_expression(src);
    if (!e) return NULL;
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    expr_free(e);
    expr_free(r);
    return s;
}

int main(int argc, char** argv) {
    const char* path = (argc > 1) ? argv[1] : "../paranoia/paranoia.m";
    char buf[2048];
    char* s;
    int pass;

    setvbuf(stdout, NULL, _IONBF, 0);
    symtab_init();
    core_init();

    if (strlen(path) > sizeof(buf) - 32) {
        printf("FAIL: path too long\n");
        return 1;
    }
    snprintf(buf, sizeof(buf), "Get[\"%s\"]; ValueQ[Paranoia`ParanoiaGate]", path);
    s = eval_str(buf);
    free(s);

    s = eval_str("Paranoia`ParanoiaGate[][\"Pass\"]");
    if (!s) {
        printf("FAIL: could not evaluate Paranoia`ParanoiaGate[] (is %s readable?)\n", path);
        return 1;
    }
    pass = (strcmp(s, "True") == 0);
    if (!pass) printf("PARANOIA gate result: %s\n", s);
    free(s);

    if (!pass) {
        printf("PARANOIA: a FAILURE or SERIOUS DEFECT was reported that is not on the known list.\n");
        return 1;
    }
    printf("All PARANOIA modes passed the gate.\n");
    return 0;
}
