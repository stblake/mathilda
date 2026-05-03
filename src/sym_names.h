/*
 * sym_names.h
 *
 * Cached, interned pointers for well-known symbol names. After
 * sym_names_init() runs, every `SYM_*` is the canonical pointer that
 * `expr_new_symbol("Foo")` and `intern_symbol("Foo")` return for that
 * name. This lets hot evaluator paths replace
 *
 *     strcmp(head->data.symbol, "List") == 0
 *
 * with the much cheaper
 *
 *     head->data.symbol == SYM_List
 *
 * Only names compared in evaluator hot paths (or otherwise pervasive)
 * need to live here; the interner makes ad-hoc pointer comparisons
 * possible everywhere else by calling intern_symbol() directly.
 */

#ifndef SYM_NAMES_H
#define SYM_NAMES_H

#ifdef __cplusplus
extern "C" {
#endif

extern const char* SYM_List;
extern const char* SYM_Plus;
extern const char* SYM_Times;
extern const char* SYM_Power;
extern const char* SYM_Sequence;
extern const char* SYM_Hold;
extern const char* SYM_HoldComplete;
extern const char* SYM_HoldPattern;
extern const char* SYM_HoldForm;
extern const char* SYM_Unevaluated;
extern const char* SYM_Evaluate;
extern const char* SYM_Set;
extern const char* SYM_SetDelayed;
extern const char* SYM_Rule;
extern const char* SYM_RuleDelayed;
extern const char* SYM_Pattern;
extern const char* SYM_Blank;
extern const char* SYM_BlankSequence;
extern const char* SYM_BlankNullSequence;
extern const char* SYM_PatternTest;
extern const char* SYM_Condition;
extern const char* SYM_Function;
extern const char* SYM_Slot;
extern const char* SYM_SlotSequence;
extern const char* SYM_CompoundExpression;
extern const char* SYM_If;
extern const char* SYM_True;
extern const char* SYM_False;
extern const char* SYM_Null;
extern const char* SYM_Rational;
extern const char* SYM_Complex;
extern const char* SYM_Part;
extern const char* SYM_Composition;
extern const char* SYM_Derivative;

/* Populate every SYM_* by interning its name string. Idempotent: safe
 * to call repeatedly. Must run before any consumer reads a SYM_*
 * pointer; in practice it is called from core_init(). */
void sym_names_init(void);

#ifdef __cplusplus
}
#endif

#endif /* SYM_NAMES_H */
