/*
 * sym_names.c
 *
 * Definitions and one-shot init for the cached SYM_* pointers declared
 * in sym_names.h.
 */

#include "sym_names.h"
#include "sym_intern.h"
#include <stddef.h>  /* NULL */

const char* SYM_List = NULL;
const char* SYM_Plus = NULL;
const char* SYM_Times = NULL;
const char* SYM_Power = NULL;
const char* SYM_Sequence = NULL;
const char* SYM_Hold = NULL;
const char* SYM_HoldComplete = NULL;
const char* SYM_HoldPattern = NULL;
const char* SYM_HoldForm = NULL;
const char* SYM_Unevaluated = NULL;
const char* SYM_Evaluate = NULL;
const char* SYM_Set = NULL;
const char* SYM_SetDelayed = NULL;
const char* SYM_Rule = NULL;
const char* SYM_RuleDelayed = NULL;
const char* SYM_Pattern = NULL;
const char* SYM_Blank = NULL;
const char* SYM_BlankSequence = NULL;
const char* SYM_BlankNullSequence = NULL;
const char* SYM_PatternTest = NULL;
const char* SYM_Condition = NULL;
const char* SYM_Function = NULL;
const char* SYM_Slot = NULL;
const char* SYM_SlotSequence = NULL;
const char* SYM_CompoundExpression = NULL;
const char* SYM_If = NULL;
const char* SYM_True = NULL;
const char* SYM_False = NULL;
const char* SYM_Null = NULL;
const char* SYM_Rational = NULL;
const char* SYM_Complex = NULL;
const char* SYM_Part = NULL;
const char* SYM_Composition = NULL;
const char* SYM_Derivative = NULL;

void sym_names_init(void) {
    /* intern_symbol is idempotent and stable, so this can run multiple
     * times without re-allocating. */
    SYM_List               = intern_symbol("List");
    SYM_Plus               = intern_symbol("Plus");
    SYM_Times              = intern_symbol("Times");
    SYM_Power              = intern_symbol("Power");
    SYM_Sequence           = intern_symbol("Sequence");
    SYM_Hold               = intern_symbol("Hold");
    SYM_HoldComplete       = intern_symbol("HoldComplete");
    SYM_HoldPattern        = intern_symbol("HoldPattern");
    SYM_HoldForm           = intern_symbol("HoldForm");
    SYM_Unevaluated        = intern_symbol("Unevaluated");
    SYM_Evaluate           = intern_symbol("Evaluate");
    SYM_Set                = intern_symbol("Set");
    SYM_SetDelayed         = intern_symbol("SetDelayed");
    SYM_Rule               = intern_symbol("Rule");
    SYM_RuleDelayed        = intern_symbol("RuleDelayed");
    SYM_Pattern            = intern_symbol("Pattern");
    SYM_Blank              = intern_symbol("Blank");
    SYM_BlankSequence      = intern_symbol("BlankSequence");
    SYM_BlankNullSequence  = intern_symbol("BlankNullSequence");
    SYM_PatternTest        = intern_symbol("PatternTest");
    SYM_Condition          = intern_symbol("Condition");
    SYM_Function           = intern_symbol("Function");
    SYM_Slot               = intern_symbol("Slot");
    SYM_SlotSequence       = intern_symbol("SlotSequence");
    SYM_CompoundExpression = intern_symbol("CompoundExpression");
    SYM_If                 = intern_symbol("If");
    SYM_True               = intern_symbol("True");
    SYM_False              = intern_symbol("False");
    SYM_Null               = intern_symbol("Null");
    SYM_Rational           = intern_symbol("Rational");
    SYM_Complex            = intern_symbol("Complex");
    SYM_Part               = intern_symbol("Part");
    SYM_Composition        = intern_symbol("Composition");
    SYM_Derivative         = intern_symbol("Derivative");
}
