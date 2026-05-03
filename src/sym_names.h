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
 * Every name compared via strcmp on a symbol field anywhere in src/
 * lives here; ad-hoc pointer comparisons elsewhere remain possible by
 * calling intern_symbol() directly.
 */

#ifndef SYM_NAMES_H
#define SYM_NAMES_H

#ifdef __cplusplus
extern "C" {
#endif

extern const char* SYM_Abort;
extern const char* SYM_Abs;
extern const char* SYM_AbsRules;
extern const char* SYM_Algebraics;
extern const char* SYM_All;
extern const char* SYM_Alternatives;
extern const char* SYM_And;
extern const char* SYM_Apart;
extern const char* SYM_ArcCos;
extern const char* SYM_ArcCosh;
extern const char* SYM_ArcCot;
extern const char* SYM_ArcCoth;
extern const char* SYM_ArcCsc;
extern const char* SYM_ArcCsch;
extern const char* SYM_ArcSec;
extern const char* SYM_ArcSech;
extern const char* SYM_ArcSin;
extern const char* SYM_ArcSinh;
extern const char* SYM_ArcTan;
extern const char* SYM_ArcTanh;
extern const char* SYM_AssumptionRules;
extern const char* SYM_Automatic;
extern const char* SYM_Base;
extern const char* SYM_BlakeRationalBaseDescent;
extern const char* SYM_Blank;
extern const char* SYM_BlankNullSequence;
extern const char* SYM_BlankSequence;
extern const char* SYM_Block;
extern const char* SYM_Boole;
extern const char* SYM_Booleans;
extern const char* SYM_Break;
extern const char* SYM_CFRAC;
extern const char* SYM_Catalan;
extern const char* SYM_Ceiling;
extern const char* SYM_CollectPerVariable;
extern const char* SYM_CompensatedSummation;
extern const char* SYM_Complex;
extern const char* SYM_ComplexInfinity;
extern const char* SYM_Complexes;
extern const char* SYM_Composites;
extern const char* SYM_Composition;
extern const char* SYM_CompoundExpression;
extern const char* SYM_Condition;
extern const char* SYM_Conjugate;
extern const char* SYM_Continue;
extern const char* SYM_Cos;
extern const char* SYM_Cosh;
extern const char* SYM_Cot;
extern const char* SYM_Coth;
extern const char* SYM_Csc;
extern const char* SYM_Csch;
extern const char* SYM_Degree;
extern const char* SYM_Derivative;
extern const char* SYM_DirectedInfinity;
extern const char* SYM_DiscreteDelta;
extern const char* SYM_Divide;
extern const char* SYM_Dixon;
extern const char* SYM_E;
extern const char* SYM_ECM;
extern const char* SYM_Equal;
extern const char* SYM_EulerGamma;
extern const char* SYM_Evaluate;
extern const char* SYM_Except;
extern const char* SYM_Exp;
extern const char* SYM_Factor;
extern const char* SYM_FactorSquareFree;
extern const char* SYM_FactorTerms;
extern const char* SYM_False;
extern const char* SYM_Fermat;
extern const char* SYM_Flat;
extern const char* SYM_Floor;
extern const char* SYM_FractionalPart;
extern const char* SYM_Frobenius;
extern const char* SYM_FromAbove;
extern const char* SYM_FromBelow;
extern const char* SYM_FullForm;
extern const char* SYM_Function;
extern const char* SYM_Gamma;
extern const char* SYM_Glaisher;
extern const char* SYM_GoldenRatio;
extern const char* SYM_Greater;
extern const char* SYM_GreaterEqual;
extern const char* SYM_Heads;
extern const char* SYM_HeavisideTheta;
extern const char* SYM_Hold;
extern const char* SYM_HoldAll;
extern const char* SYM_HoldAllComplete;
extern const char* SYM_HoldComplete;
extern const char* SYM_HoldFirst;
extern const char* SYM_HoldForm;
extern const char* SYM_HoldPattern;
extern const char* SYM_HoldRest;
extern const char* SYM_I;
extern const char* SYM_Identity;
extern const char* SYM_If;
extern const char* SYM_Im;
extern const char* SYM_Implies;
extern const char* SYM_Indeterminate;
extern const char* SYM_Infinity;
extern const char* SYM_InputForm;
extern const char* SYM_Integer;
extern const char* SYM_IntegerPart;
extern const char* SYM_Integers;
extern const char* SYM_Inverse;
extern const char* SYM_InverseFunction;
extern const char* SYM_Khinchin;
extern const char* SYM_KroneckerDelta;
extern const char* SYM_Less;
extern const char* SYM_LessEqual;
extern const char* SYM_List;
extern const char* SYM_Listable;
extern const char* SYM_Locked;
extern const char* SYM_Log;
extern const char* SYM_Log1p;
extern const char* SYM_LogExpRules;
extern const char* SYM_Longest;
extern const char* SYM_MachineEpsilon;
extern const char* SYM_MachinePrecision;
extern const char* SYM_Method;
extern const char* SYM_Mod;
extern const char* SYM_Module;
extern const char* SYM_Modulus;
extern const char* SYM_NHoldRest;
extern const char* SYM_None;
extern const char* SYM_Not;
extern const char* SYM_Null;
extern const char* SYM_NumericFunction;
extern const char* SYM_OneIdentity;
extern const char* SYM_Optional;
extern const char* SYM_OptionsPattern;
extern const char* SYM_Or;
extern const char* SYM_Orderless;
extern const char* SYM_Overflow;
extern const char* SYM_Part;
extern const char* SYM_Pattern;
extern const char* SYM_PatternTest;
extern const char* SYM_Pi;
extern const char* SYM_Piecewise;
extern const char* SYM_Plus;
extern const char* SYM_PollardRho;
extern const char* SYM_Power;
extern const char* SYM_Primes;
extern const char* SYM_Protected;
extern const char* SYM_Quit;
extern const char* SYM_Quotient;
extern const char* SYM_Rational;
extern const char* SYM_Rationals;
extern const char* SYM_Re;
extern const char* SYM_ReadProtected;
extern const char* SYM_Real;
extern const char* SYM_Reals;
extern const char* SYM_Repeated;
extern const char* SYM_RepeatedNull;
extern const char* SYM_Return;
extern const char* SYM_Round;
extern const char* SYM_Rule;
extern const char* SYM_RuleDelayed;
extern const char* SYM_SameQ;
extern const char* SYM_SameTest;
extern const char* SYM_Sec;
extern const char* SYM_Sech;
extern const char* SYM_Sequence;
extern const char* SYM_SequenceHold;
extern const char* SYM_SeriesData;
extern const char* SYM_Set;
extern const char* SYM_SetDelayed;
extern const char* SYM_ShanksSquareForms;
extern const char* SYM_Shortest;
extern const char* SYM_Sign;
extern const char* SYM_Sin;
extern const char* SYM_Sinh;
extern const char* SYM_Slot;
extern const char* SYM_SlotSequence;
extern const char* SYM_Span;
extern const char* SYM_Sqrt;
extern const char* SYM_String;
extern const char* SYM_Symbol;
extern const char* SYM_Table;
extern const char* SYM_Tan;
extern const char* SYM_Tanh;
extern const char* SYM_TeXForm;
extern const char* SYM_Temporary;
extern const char* SYM_Throw;
extern const char* SYM_Times;
extern const char* SYM_TrialDivision;
extern const char* SYM_TrigExpand;
extern const char* SYM_TrigFactor;
extern const char* SYM_TrigReduce;
extern const char* SYM_TrigRoundtrip;
extern const char* SYM_TrigToExp;
extern const char* SYM_True;
extern const char* SYM_TwoSided;
extern const char* SYM_Unequal;
extern const char* SYM_Unevaluated;
extern const char* SYM_UnitStep;
extern const char* SYM_UnsameQ;
extern const char* SYM_UpTo;
extern const char* SYM_Verbatim;
extern const char* SYM_With;
extern const char* SYM_Xor;

/* Populate every SYM_* by interning its name string. Idempotent: safe
 * to call repeatedly. Must run before any consumer reads a SYM_*
 * pointer; in practice it is called from core_init(). */
void sym_names_init(void);

#ifdef __cplusplus
}
#endif

#endif /* SYM_NAMES_H */
