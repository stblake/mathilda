#include "attr.h"
#include "symtab.h"
#include "eval.h"
#include "sym_names.h"
#include "message.h"   /* mth_message */
#include "print.h"     /* expr_to_string */
#include "common.h"    /* head_is */
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

// Simple attribute dictionary for built-ins
typedef struct {
    const char* name;
    uint32_t attrs;
} SymbolAttr;

static SymbolAttr builtin_attrs[] = {
    {"Hold", ATTR_HOLDALL | ATTR_PROTECTED},
    {"HoldFirst", ATTR_HOLDFIRST | ATTR_PROTECTED},
    {"HoldRest", ATTR_HOLDREST | ATTR_PROTECTED},
    {"HoldComplete", ATTR_HOLDALLCOMPLETE | ATTR_PROTECTED},
    {"HoldPattern", ATTR_HOLDALL | ATTR_PROTECTED},
    {"Literal", ATTR_HOLDALL | ATTR_PROTECTED},
    {"Unevaluated", ATTR_HOLDALLCOMPLETE | ATTR_PROTECTED},
    {"Set", ATTR_HOLDFIRST | ATTR_PROTECTED | ATTR_SEQUENCEHOLD},
    {"SetDelayed", ATTR_HOLDALL | ATTR_PROTECTED | ATTR_SEQUENCEHOLD},
    {"MessageName", ATTR_HOLDFIRST | ATTR_PROTECTED},
    {"Clear", ATTR_HOLDALL | ATTR_PROTECTED},
    {"AppendTo", ATTR_HOLDFIRST | ATTR_PROTECTED},
    {"PrependTo", ATTR_HOLDFIRST | ATTR_PROTECTED},
    {"Plus", ATTR_FLAT | ATTR_ORDERLESS | ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_ONEIDENTITY | ATTR_PROTECTED},
    {"Times", ATTR_FLAT | ATTR_ORDERLESS | ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_ONEIDENTITY | ATTR_PROTECTED},
    {"Divide", ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED},
    {"Subtract", ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_ONEIDENTITY | ATTR_PROTECTED},
    {"Power", ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_ONEIDENTITY | ATTR_PROTECTED},
    {"Sqrt", ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED},
    {"Rational", ATTR_PROTECTED},
    {"Attributes", ATTR_HOLDALL | ATTR_LISTABLE | ATTR_PROTECTED},
    {"SetAttributes", ATTR_HOLDFIRST | ATTR_PROTECTED},
    {"ClearAttributes", ATTR_HOLDFIRST | ATTR_PROTECTED},
    {"OwnValues", ATTR_HOLDALL | ATTR_PROTECTED},
    {"DownValues", ATTR_HOLDALL | ATTR_PROTECTED},
    {"FormatValues", ATTR_HOLDALL | ATTR_PROTECTED},
    {"NValues", ATTR_HOLDALL | ATTR_PROTECTED},
    {"DefaultValues", ATTR_HOLDALL | ATTR_PROTECTED},
    {"Default", ATTR_PROTECTED},
    {"Format", ATTR_PROTECTED},
    {"MakeBoxes", ATTR_HOLDALLCOMPLETE},
#if UP_VALUES
    {"UpValues", ATTR_HOLDALL | ATTR_PROTECTED},
    {"UpSet", ATTR_HOLDFIRST | ATTR_PROTECTED | ATTR_SEQUENCEHOLD},
    {"UpSetDelayed", ATTR_HOLDALL | ATTR_PROTECTED | ATTR_SEQUENCEHOLD},
    {"TagSet", ATTR_HOLDALL | ATTR_PROTECTED | ATTR_SEQUENCEHOLD},
    {"TagSetDelayed", ATTR_HOLDALL | ATTR_PROTECTED | ATTR_SEQUENCEHOLD},
    {"TagUnset", ATTR_HOLDALL | ATTR_PROTECTED | ATTR_SEQUENCEHOLD},
    {"Definition", ATTR_HOLDALL | ATTR_PROTECTED},
#endif
    {"Out", ATTR_LISTABLE | ATTR_NHOLDFIRST | ATTR_PROTECTED},
    {"Overflow", ATTR_PROTECTED},
    {"Underflow", ATTR_PROTECTED},
    {"Apply", ATTR_PROTECTED},
    {"Map", ATTR_PROTECTED},
    {"MapAll", ATTR_PROTECTED},
    {"Select", ATTR_PROTECTED},
    {"FreeQ", ATTR_PROTECTED},
    {"Sort", ATTR_PROTECTED},
    {"PolynomialQ", ATTR_PROTECTED},
    {"Variables", ATTR_PROTECTED},
    {"Level", ATTR_PROTECTED},
    {"Depth", ATTR_PROTECTED},
    {"Function", ATTR_HOLDALL | ATTR_PROTECTED},
    {"MatchQ", ATTR_PROTECTED},
    {"Pattern", ATTR_HOLDFIRST | ATTR_PROTECTED},
    {"PatternTest", ATTR_HOLDREST | ATTR_PROTECTED},
    {"PatternSequence", ATTR_PROTECTED},
    {"Repeated", ATTR_PROTECTED},
    {"RepeatedNull", ATTR_PROTECTED},
    {"OrderlessPatternSequence", ATTR_PROTECTED},
    {"Verbatim", ATTR_PROTECTED},
    {"Condition", ATTR_HOLDREST | ATTR_PROTECTED},
    {"Rule", ATTR_PROTECTED | ATTR_SEQUENCEHOLD},
    {"RuleDelayed", ATTR_HOLDREST | ATTR_PROTECTED | ATTR_SEQUENCEHOLD},
    {"Equal", ATTR_PROTECTED},
    {"Unequal", ATTR_PROTECTED},
    {"SameQ", ATTR_PROTECTED},
    {"UnsameQ", ATTR_PROTECTED},
    {"Less", ATTR_PROTECTED},
    {"Greater", ATTR_PROTECTED},
    {"LessEqual", ATTR_PROTECTED},
    {"GreaterEqual", ATTR_PROTECTED},
    {"And", ATTR_FLAT | ATTR_HOLDALL | ATTR_ONEIDENTITY | ATTR_PROTECTED},
    {"Or", ATTR_FLAT | ATTR_HOLDALL | ATTR_ONEIDENTITY | ATTR_PROTECTED},
    {"Not", ATTR_PROTECTED},
    {"Xor", ATTR_FLAT | ATTR_ORDERLESS | ATTR_ONEIDENTITY | ATTR_PROTECTED},
    /* Not Flat: Equivalent[a, Equivalent[b, c]] tests a against the truth
     * value of b <=> c, which is not Equivalent[a, b, c] (audit Equivalent#1). */
    {"Equivalent", ATTR_ORDERLESS | ATTR_PROTECTED},
    {"Implies", ATTR_PROTECTED},
    {"Boole", ATTR_LISTABLE | ATTR_PROTECTED},
    {"ConditionalExpression", ATTR_PROTECTED},
    {"CompoundExpression", ATTR_HOLDALL | ATTR_PROTECTED},
    {"Table", ATTR_HOLDALL | ATTR_PROTECTED},
    {"Module", ATTR_HOLDALL | ATTR_PROTECTED},
    {"Block", ATTR_HOLDALL | ATTR_PROTECTED},
    {"With", ATTR_HOLDALL | ATTR_PROTECTED},
    {"Unique", ATTR_PROTECTED},
    {"Range", ATTR_LISTABLE | ATTR_PROTECTED},
    {"Return", ATTR_PROTECTED},
    {"Array", ATTR_PROTECTED},
    {"ConstantArray", ATTR_PROTECTED},
    {"ArrayFlatten", ATTR_PROTECTED},
    {"ArrayReshape", ATTR_PROTECTED},
    {"ArrayPad", ATTR_PROTECTED},
    {"Take", ATTR_NHOLDREST | ATTR_PROTECTED},
    {"Drop", ATTR_NHOLDREST | ATTR_PROTECTED},
    {"Extract", ATTR_NHOLDREST | ATTR_PROTECTED},
    {"Flatten", ATTR_PROTECTED},
    {"FlattenAt", ATTR_PROTECTED},
    {"Partition", ATTR_PROTECTED},
    {"RotateLeft", ATTR_PROTECTED},
    {"RotateRight", ATTR_PROTECTED},
    {"Reverse", ATTR_PROTECTED},
    {"Transpose", ATTR_PROTECTED},
    {"Mod", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Quotient", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"QuotientRemainder", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"GCD", ATTR_PROTECTED | ATTR_LISTABLE | ATTR_FLAT | ATTR_ORDERLESS | ATTR_ONEIDENTITY},
    {"LCM", ATTR_PROTECTED | ATTR_LISTABLE | ATTR_FLAT | ATTR_ORDERLESS | ATTR_ONEIDENTITY},
    {"ExtendedGCD", ATTR_PROTECTED | ATTR_LISTABLE},
    {"PrimeQ", ATTR_PROTECTED | ATTR_LISTABLE},
    {"PrimePowerQ", ATTR_PROTECTED | ATTR_LISTABLE},
    {"PossibleZeroQ", ATTR_PROTECTED},
    {"PrimePi", ATTR_PROTECTED | ATTR_LISTABLE},
    {"FactorInteger", ATTR_PROTECTED | ATTR_LISTABLE},
    {"NextPrime", ATTR_PROTECTED},        /* not Listable */
    {"Fibonacci", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"LucasL", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Re", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Im", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"ReIm", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Abs", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Conjugate", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Arg", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Sin", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Cos", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Tan", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Cot", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Sec", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Csc", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"ArcSin", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"ArcCos", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"ArcTan", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"ArcCot", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"ArcSec", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"ArcCsc", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Log", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Log10", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Log2", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Exp", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"E", ATTR_PROTECTED},
    {"I", ATTR_PROTECTED},
    {"Floor", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Ceiling", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Round", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"IntegerPart", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"IntegerDigits", ATTR_PROTECTED | ATTR_LISTABLE},
    {"IntegerLength", ATTR_PROTECTED | ATTR_LISTABLE},
    {"BitLength", ATTR_PROTECTED | ATTR_LISTABLE},
    {"BitXor", ATTR_PROTECTED | ATTR_LISTABLE | ATTR_FLAT | ATTR_ORDERLESS | ATTR_ONEIDENTITY},
    {"BitAnd", ATTR_PROTECTED | ATTR_LISTABLE | ATTR_FLAT | ATTR_ORDERLESS | ATTR_ONEIDENTITY},
    {"BitOr", ATTR_PROTECTED | ATTR_LISTABLE | ATTR_FLAT | ATTR_ORDERLESS | ATTR_ONEIDENTITY},
    {"BitNot", ATTR_PROTECTED | ATTR_LISTABLE},
    {"BitSet", ATTR_PROTECTED | ATTR_LISTABLE},
    {"BitClear", ATTR_PROTECTED | ATTR_LISTABLE},
    {"BitGet", ATTR_PROTECTED | ATTR_LISTABLE},
    {"BitShiftLeft", ATTR_PROTECTED | ATTR_LISTABLE},
    {"BitShiftRight", ATTR_PROTECTED | ATTR_LISTABLE},
    {"DigitCount", ATTR_PROTECTED},
    {"FromDigits", ATTR_PROTECTED},
    {"FractionalPart", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Timing", ATTR_HOLDALL | ATTR_PROTECTED | ATTR_SEQUENCEHOLD},
    {"AbsoluteTiming", ATTR_HOLDALL | ATTR_PROTECTED | ATTR_SEQUENCEHOLD},
    {"EchoTiming", ATTR_HOLDFIRST | ATTR_PROTECTED | ATTR_SEQUENCEHOLD},
    {"Echo", ATTR_PROTECTED},
    {"RepeatedTiming", ATTR_HOLDFIRST | ATTR_PROTECTED | ATTR_SEQUENCEHOLD},
    {"TimeConstrained", ATTR_HOLDALL | ATTR_PROTECTED},
    {"MemoryConstrained", ATTR_HOLDALL | ATTR_PROTECTED},
    {"Dot", ATTR_FLAT | ATTR_ONEIDENTITY | ATTR_PROTECTED},
    {"NonCommutativeMultiply", ATTR_FLAT | ATTR_ONEIDENTITY | ATTR_PROTECTED},
    {"Det", ATTR_PROTECTED},
    {"Cross", ATTR_PROTECTED},
    {"Norm", ATTR_PROTECTED},
    {"Sequence", ATTR_PROTECTED},
    {"Nothing", ATTR_PROTECTED},
    {"SequenceHold", ATTR_PROTECTED},

    /* System symbols that carried no attributes at all (builtin audit, S7:
     * List#1, Blank#1, Head#1, Append#1, First#1, ...).  Every System symbol
     * that is not a user-settable variable is Protected, as in the reference:
     * an unprotected `List = 5` turned every list in the session into
     * 5[...].  (Locked is deliberately not modelled; see attr.h.)  Settable
     * on purpose and so absent here: the $-variables ($Context, $Path,
     * $RecursionLimit, ...), MakeBoxes and Derivative (users attach rules to
     * them), and Default (assignments to it go through Set's own check). */
    {"List", ATTR_PROTECTED},
    {"Join", ATTR_FLAT | ATTR_ONEIDENTITY | ATTR_PROTECTED},
    {"Head", ATTR_PROTECTED},
    {"Length", ATTR_PROTECTED},
    {"Dimensions", ATTR_PROTECTED},
    {"First", ATTR_PROTECTED},
    {"Last", ATTR_PROTECTED},
    {"Most", ATTR_PROTECTED},
    {"Rest", ATTR_PROTECTED},
    {"Append", ATTR_PROTECTED},
    {"Prepend", ATTR_PROTECTED},
    {"Insert", ATTR_PROTECTED},
    {"Delete", ATTR_PROTECTED},
    {"Through", ATTR_PROTECTED},
    {"Blank", ATTR_PROTECTED},
    {"BlankSequence", ATTR_PROTECTED},
    {"BlankNullSequence", ATTR_PROTECTED},
    {"Alternatives", ATTR_PROTECTED},
    {"Optional", ATTR_PROTECTED},
    {"OptionsPattern", ATTR_PROTECTED},
    {"Longest", ATTR_PROTECTED},
    {"Shortest", ATTR_PROTECTED},
    {"UpTo", ATTR_PROTECTED},
    {"Span", ATTR_PROTECTED},
    {"Integer", ATTR_PROTECTED},
    {"Real", ATTR_PROTECTED},
    {"Complex", ATTR_PROTECTED},
    {"String", ATTR_PROTECTED},
    {"True", ATTR_PROTECTED},
    {"False", ATTR_PROTECTED},
    {"Null", ATTR_PROTECTED},
    {"Infinity", ATTR_PROTECTED},
    {"All", ATTR_PROTECTED},
    {"None", ATTR_PROTECTED},
    {"Automatic", ATTR_PROTECTED},
    {"$Failed", ATTR_PROTECTED},
    {"Integers", ATTR_PROTECTED},
    {"Rationals", ATTR_PROTECTED},
    {"Reals", ATTR_PROTECTED},
    {"Complexes", ATTR_PROTECTED},
    {"Flat", ATTR_PROTECTED},
    {"Listable", ATTR_PROTECTED},
    {"HoldAll", ATTR_PROTECTED},
    {"HoldAllComplete", ATTR_PROTECTED},
    {"OneIdentity", ATTR_PROTECTED},
    {"Orderless", ATTR_PROTECTED},
    {"Protected", ATTR_PROTECTED},
    {"Temporary", ATTR_PROTECTED},
    {"Number", ATTR_PROTECTED},
    {"Character", ATTR_PROTECTED},
    {"Word", ATTR_PROTECTED},
    {"Record", ATTR_PROTECTED},
    {"Expression", ATTR_PROTECTED},
    {"Byte", ATTR_PROTECTED},
    {"EndOfFile", ATTR_PROTECTED},
    {"Cubics", ATTR_PROTECTED},
    {"Quartics", ATTR_PROTECTED},
    {"ChessboardDistance", ATTR_PROTECTED},
    /* Option names. */
    {"AccuracyGoal", ATTR_PROTECTED},
    {"PrecisionGoal", ATTR_PROTECTED},
    {"AspectRatio", ATTR_PROTECTED},
    {"Assumptions", ATTR_PROTECTED},
    {"Frame", ATTR_PROTECTED},
    {"FrameStyle", ATTR_PROTECTED},
    {"FrameTicks", ATTR_PROTECTED},
    {"ImageSize", ATTR_PROTECTED},
    {"Method", ATTR_PROTECTED},
    {"Modulus", ATTR_PROTECTED},
    {"TimeConstraint", ATTR_PROTECTED},
    {"VerifySolutions", ATTR_PROTECTED},
    {"GeneratedParameters", ATTR_PROTECTED},
    {"InverseFunctions", ATTR_PROTECTED},
    {"TransformationFunctions", ATTR_PROTECTED},
    {"NullRecords", ATTR_PROTECTED},
    {"NullWords", ATTR_PROTECTED},
    {"RecordSeparators", ATTR_PROTECTED},
    {"TokenWords", ATTR_PROTECTED},
    {"WordSeparators", ATTR_PROTECTED},
    {NULL, 0}
};

/* Effective attributes of a resolved SymbolDef: an O(1) field read, no
 * symbol-table lookup. The evaluator resolves the head's def once per step and
 * calls this. builtin_attrs[] above is only the initial seed (folded into
 * def->attributes by attr_init); after that the field is the single source of
 * truth, so Unprotect / ClearAttributes on a builtin really do take effect. */
uint32_t get_attributes_def(SymbolDef* def) {
    return def ? def->attributes : ATTR_NONE;
}

/* Current attributes of a call's symbol head, through the def the symbol node
 * caches (the evaluator fills it on its first pass over the node, so for a
 * builtin it is a couple of loads; a node built in C pays one lookup, once).
 * For a builtin that must honour a user's SetAttributes / ClearAttributes on
 * its own head -- Plus and Times re-check Orderless before sorting. */
uint32_t call_head_attributes(Expr* call) {
    if (!call || call->type != EXPR_FUNCTION) return ATTR_NONE;
    Expr* h = call->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return ATTR_NONE;
    SymbolDef* d = h->data.symbol.def;
    if (!d) {
        d = symtab_get_def(h->data.symbol.name);
        h->data.symbol.def = d;      /* benign metadata, as in evaluate_step */
    }
    return get_attributes_def(d);
}

uint32_t get_attributes(const char* symbol_name) {
    if (!symbol_name) return ATTR_NONE;
    return get_attributes_def(symtab_get_def(symbol_name));
}

void set_attributes(const char* symbol_name, uint32_t attrs) {
    if (!symbol_name) return;
    SymbolDef* def = symtab_get_def(symbol_name);
    if (def) {
        if (def->attributes != attrs) {
            def->attributes = attrs;
            /* Attribute changes (Hold*, Listable, Flat, Orderless,
             * NumericFunction, Protected, OneIdentity, ...) all change
             * how the evaluator handles a head, so a cached evaluation
             * could be stale. Advance the rule epoch (and eval clock) to
             * invalidate -- including any GROUND node using this head. */
            eval_rule_epoch_bump();
        }
    }
}

static uint32_t string_to_attribute(const char* name) {
    if (strcmp(name, "HoldFirst") == 0) return ATTR_HOLDFIRST;
    if (strcmp(name, "HoldRest") == 0) return ATTR_HOLDREST;
    if (strcmp(name, "HoldAll") == 0) return ATTR_HOLDALL;
    if (strcmp(name, "HoldAllComplete") == 0) return ATTR_HOLDALLCOMPLETE;
    if (strcmp(name, "Flat") == 0) return ATTR_FLAT;
    if (strcmp(name, "Orderless") == 0) return ATTR_ORDERLESS;
    if (strcmp(name, "Listable") == 0) return ATTR_LISTABLE;
    if (strcmp(name, "NumericFunction") == 0) return ATTR_NUMERICFUNCTION;
    if (strcmp(name, "Protected") == 0) return ATTR_PROTECTED;
    if (strcmp(name, "OneIdentity") == 0) return ATTR_ONEIDENTITY;
    if (strcmp(name, "NHoldAll") == 0) return ATTR_NHOLDALL;
    if (strcmp(name, "NHoldFirst") == 0) return ATTR_NHOLDFIRST;
    if (strcmp(name, "NHoldRest") == 0) return ATTR_NHOLDREST;
    if (strcmp(name, "Temporary") == 0) return ATTR_TEMPORARY;
    if (strcmp(name, "SequenceHold") == 0) return ATTR_SEQUENCEHOLD;
    if (strcmp(name, "Constant") == 0) return ATTR_CONSTANT;
    return ATTR_NONE;
}

/* Attributes[s] = spec: replace s's attributes by exactly those named in
 * spec (one name or a List of names).  False, changing nothing, when spec
 * names something that is not an attribute. */
bool attr_assign(const char* symbol_name, Expr* spec) {
    if (!symbol_name || !spec) return false;
    Expr* one[1] = { spec };
    Expr** items = one;
    size_t n = 1;
    if (spec->type == EXPR_FUNCTION && spec->data.function.head->type == EXPR_SYMBOL
        && spec->data.function.head->data.symbol.name == SYM_List) {
        items = spec->data.function.args;
        n = spec->data.function.arg_count;
    }
    uint32_t bits = ATTR_NONE;
    for (size_t i = 0; i < n; i++) {
        const char* nm = items[i]->type == EXPR_SYMBOL ? items[i]->data.symbol.name
                       : items[i]->type == EXPR_STRING ? items[i]->data.string : NULL;
        uint32_t b = nm ? string_to_attribute(nm) : ATTR_NONE;
        if (b == ATTR_NONE) return false;
        bits |= b;
    }
    set_attributes(symbol_name, bits);
    return true;
}

static void remove_single_attribute(SymbolDef* def, Expr* attr_expr) {
    const char* attr_name = NULL;
    if (attr_expr->type == EXPR_SYMBOL) attr_name = attr_expr->data.symbol.name;
    else if (attr_expr->type == EXPR_STRING) attr_name = attr_expr->data.string;

    if (attr_name) {
        uint32_t bit = string_to_attribute(attr_name);
        /* HoldAll (NHoldAll) is one attribute stored as the two bits of
         * HoldFirst | HoldRest.  Clearing the composite needs the whole
         * composite present (ClearAttributes[h, HoldAll] on a HoldFirst h
         * keeps HoldFirst), and clearing one half of a full composite is a
         * no-op (ClearAttributes[g, HoldFirst] on a HoldAll g keeps HoldAll)
         * -- audit ClearAttributes#1/#2. */
        static const uint32_t composites[2] = { ATTR_HOLDALL, ATTR_NHOLDALL };
        for (int k = 0; k < 2; k++) {
            uint32_t c = composites[k];
            if (bit == c && (def->attributes & c) != c) return;
            if (bit != c && (bit & c) == bit && (def->attributes & c) == c) return;
        }
        if (bit != ATTR_NONE && (def->attributes & bit)) {
            def->attributes &= ~bit;
            /* Real attribute change -- invalidate cached evaluations (and any
             * GROUND node using this head) by advancing the rule epoch. */
            eval_rule_epoch_bump();
        }
    }
}

static void add_single_attribute(SymbolDef* def, Expr* attr_expr) {
    const char* attr_name = NULL;
    if (attr_expr->type == EXPR_SYMBOL) attr_name = attr_expr->data.symbol.name;
    else if (attr_expr->type == EXPR_STRING) attr_name = attr_expr->data.string;

    if (attr_name) {
        uint32_t bit = string_to_attribute(attr_name);
        /* A composite (HoldAll = HoldFirst | HoldRest) is added unless ALL of
         * it is present: SetAttributes[h, HoldAll] on a HoldFirst h used to
         * see the overlap and add nothing. */
        if (bit != ATTR_NONE && (def->attributes & bit) != bit) {
            def->attributes |= bit;
            /* Real attribute change -- invalidate cached evaluations (and any
             * GROUND node using this head) by advancing the rule epoch. */
            eval_rule_epoch_bump();
        }
    }
}

/* Is `name` an attribute the language knows but Mathilda does not store?
 * ReadProtected (deliberately absent: ?f always shows definitions), Locked and
 * Stub.  SetAttributes / ClearAttributes accept them as no-ops rather than
 * reporting them as unknown (packages routinely set ReadProtected). */
static bool attr_known_unstored(const char* name) {
    return strcmp(name, "ReadProtected") == 0 || strcmp(name, "Locked") == 0
        || strcmp(name, "Stub") == 0;
}

/* Validate an attribute spec (one name or a List of names) for head `head`
 * (SetAttributes / ClearAttributes).  An element that is not a known attribute
 * emits head::attnf ("Foo is not a known attribute.") and the whole call is
 * refused, before anything changes (audit SetAttributes#1, ClearAttributes#3). */
static bool attr_spec_valid(const char* head, Expr* attr_spec) {
    Expr* one[1] = { attr_spec };
    Expr** items = one;
    size_t n = 1;
    if (head_is(attr_spec, SYM_List)) {
        items = attr_spec->data.function.args;
        n = attr_spec->data.function.arg_count;
    }
    for (size_t i = 0; i < n; i++) {
        const char* nm = items[i]->type == EXPR_SYMBOL ? items[i]->data.symbol.name
                       : items[i]->type == EXPR_STRING ? items[i]->data.string : NULL;
        if (nm && (string_to_attribute(nm) != ATTR_NONE || attr_known_unstored(nm)))
            continue;
        char* es = expr_to_string(items[i]);
        mth_message(head, "attnf", "%s is not a known attribute.", es ? es : "");
        free(es);
        return false;
    }
    return true;
}

/* The symbol name a SetAttributes / ClearAttributes target names: a symbol or
 * a string.  Anything else emits head::sym ("Argument e at position 1 is
 * expected to be a symbol.") and gives NULL (audit SetAttributes#2,
 * ClearAttributes#4). */
static const char* attr_target_name(const char* head, Expr* sym_expr) {
    if (sym_expr->type == EXPR_SYMBOL) return sym_expr->data.symbol.name;
    if (sym_expr->type == EXPR_STRING) return sym_expr->data.string;
    char* es = expr_to_string(sym_expr);
    mth_message(head, "sym", "Argument %s at position 1 is expected to be a symbol.",
                es ? es : "");
    free(es);
    return NULL;
}

static void set_attributes_for_symbol(const char* sym_name, Expr* attr_spec) {
    SymbolDef* def = symtab_get_def(sym_name);
    if (!def) return;

    // Attribute spec can be a single attribute or a list of attributes
    if (attr_spec->type == EXPR_SYMBOL || attr_spec->type == EXPR_STRING) {
        add_single_attribute(def, attr_spec);
    } else if (head_is(attr_spec, SYM_List)) {
        for (size_t i = 0; i < attr_spec->data.function.arg_count; i++) {
            add_single_attribute(def, attr_spec->data.function.args[i]);
        }
    }
}

static void clear_attributes_for_symbol(const char* sym_name, Expr* attr_spec) {
    /* Never create the symbol: clearing an attribute of a symbol that does
     * not exist has nothing to clear (audit ClearAttributes#6). */
    SymbolDef* def = symtab_lookup(sym_name);
    if (!def) return;

    // Attribute spec can be a single attribute or a list of attributes
    if (attr_spec->type == EXPR_SYMBOL || attr_spec->type == EXPR_STRING) {
        remove_single_attribute(def, attr_spec);
    } else if (head_is(attr_spec, SYM_List)) {
        for (size_t i = 0; i < attr_spec->data.function.arg_count; i++) {
            remove_single_attribute(def, attr_spec->data.function.args[i]);
        }
    }
}

/* Shared driver of SetAttributes / ClearAttributes: the target is a symbol, a
 * string, or a List of them; the spec is validated first.  An invalid target
 * or attribute leaves the call unevaluated. */
static Expr* attr_set_or_clear(Expr* res, bool setting) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    const char* head = setting ? "SetAttributes" : "ClearAttributes";
    Expr* sym_spec = res->data.function.args[0];
    Expr* attr_spec = res->data.function.args[1];

    Expr* one[1] = { sym_spec };
    Expr** targets = one;
    size_t nt = 1;
    if (head_is(sym_spec, SYM_List)) {
        targets = sym_spec->data.function.args;
        nt = sym_spec->data.function.arg_count;
    }
    for (size_t i = 0; i < nt; i++)
        if (!attr_target_name(head, targets[i])) return NULL;
    if (!attr_spec_valid(head, attr_spec)) return NULL;

    for (size_t i = 0; i < nt; i++) {
        const char* nm = targets[i]->type == EXPR_SYMBOL ? targets[i]->data.symbol.name
                                                         : targets[i]->data.string;
        if (setting) set_attributes_for_symbol(nm, attr_spec);
        else         clear_attributes_for_symbol(nm, attr_spec);
    }
    return expr_new_symbol(SYM_Null);
}

Expr* builtin_clear_attributes(Expr* res) {
    return attr_set_or_clear(res, false);
}

Expr* builtin_set_attributes(Expr* res) {
    return attr_set_or_clear(res, true);
}

static Expr* attributes_of(Expr* arg);

/* Attributes::ssle for an argument that names no symbol. */
static void attributes_ssle(const Expr* el) {
    char* es = expr_to_string((Expr*)el);
    mth_message("Attributes", "ssle",
                "Symbol, string, or HoldPattern[symbol] expected at position 1 in Attributes[%s].",
                es ? es : "");
    free(es);
}

/* Attributes[s] / Attributes["name"] -> the attribute list of one symbol;
 * Attributes[{s1, s2, ...}] -> {Attributes[s1], Attributes[s2], ...} (the
 * DifferentialForms package tests FreeQ[Attributes[{x}], Constant]). */
Expr* builtin_attributes(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol.name == SYM_List) {
        size_t n = arg->data.function.arg_count;
        Expr** items = malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t k = 0; k < n; k++) {
            Expr* el = arg->data.function.args[k];
            items[k] = attributes_of(el);
            if (!items[k]) {
                /* The bad element stays Attributes[el] (with ::ssle) and
                 * the rest are still answered. */
                attributes_ssle(el);
                Expr* one = expr_copy(el);
                items[k] = expr_new_function(expr_new_symbol("Attributes"), &one, 1);
            }
        }
        Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, n);
        free(items);
        return out;
    }
    Expr* out = attributes_of(arg);
    if (!out) attributes_ssle(arg);   /* audit Attributes#2 */
    return out;
}

static Expr* attributes_of(Expr* arg) {
    const char* sym_name = NULL;
    if (head_is(arg, SYM_HoldPattern) && arg->data.function.arg_count == 1)
        arg = arg->data.function.args[0];   /* Attributes[HoldPattern[s]] */
    if (arg->type == EXPR_SYMBOL) sym_name = arg->data.symbol.name;
    else if (arg->type == EXPR_STRING) sym_name = arg->data.string;
    
    if (!sym_name) return NULL;
    
    uint32_t attrs = get_attributes(sym_name);
    
    // Count attributes
    size_t count = 0;
    if (attrs & ATTR_CONSTANT) count++;
    if (attrs & ATTR_HOLDFIRST) count++;
    if (attrs & ATTR_HOLDREST) count++;
    if (attrs & ATTR_HOLDALLCOMPLETE) count++;
    if (attrs & ATTR_FLAT) count++;
    if (attrs & ATTR_ORDERLESS) count++;
    if (attrs & ATTR_LISTABLE) count++;
    if (attrs & ATTR_NUMERICFUNCTION) count++;
    if (attrs & ATTR_PROTECTED) count++;
    if (attrs & ATTR_ONEIDENTITY) count++;
    if ((attrs & ATTR_NHOLDALL) == ATTR_NHOLDALL) {
        count++;
    } else {
        if (attrs & ATTR_NHOLDFIRST) count++;
        if (attrs & ATTR_NHOLDREST) count++;
    }
    if (attrs & ATTR_TEMPORARY) count++;
    if (attrs & ATTR_SEQUENCEHOLD) count++;

    Expr** attr_list = malloc(sizeof(Expr*) * count);
    size_t i = 0;
    if (attrs & ATTR_CONSTANT) attr_list[i++] = expr_new_symbol(SYM_Constant);
    if (attrs & ATTR_FLAT) attr_list[i++] = expr_new_symbol(SYM_Flat);
    if ((attrs & ATTR_HOLDALL) == ATTR_HOLDALL) {
        attr_list[i++] = expr_new_symbol(SYM_HoldAll);
    } else {
        if (attrs & ATTR_HOLDFIRST) attr_list[i++] = expr_new_symbol(SYM_HoldFirst);
        if (attrs & ATTR_HOLDREST) attr_list[i++] = expr_new_symbol(SYM_HoldRest);
    }
    if (attrs & ATTR_HOLDALLCOMPLETE) attr_list[i++] = expr_new_symbol(SYM_HoldAllComplete);
    if (attrs & ATTR_LISTABLE) attr_list[i++] = expr_new_symbol(SYM_Listable);
    /* Alphabetical, case-insensitively: NHold* sorts before NumericFunction. */
    if ((attrs & ATTR_NHOLDALL) == ATTR_NHOLDALL) {
        attr_list[i++] = expr_new_symbol(SYM_NHoldAll);
    } else {
        if (attrs & ATTR_NHOLDFIRST) attr_list[i++] = expr_new_symbol(SYM_NHoldFirst);
        if (attrs & ATTR_NHOLDREST) attr_list[i++] = expr_new_symbol(SYM_NHoldRest);
    }
    if (attrs & ATTR_NUMERICFUNCTION) attr_list[i++] = expr_new_symbol(SYM_NumericFunction);
    if (attrs & ATTR_ONEIDENTITY) attr_list[i++] = expr_new_symbol(SYM_OneIdentity);
    if (attrs & ATTR_ORDERLESS) attr_list[i++] = expr_new_symbol(SYM_Orderless);
    if (attrs & ATTR_PROTECTED) attr_list[i++] = expr_new_symbol(SYM_Protected);
    if (attrs & ATTR_SEQUENCEHOLD) attr_list[i++] = expr_new_symbol(SYM_SequenceHold);
    if (attrs & ATTR_TEMPORARY) attr_list[i++] = expr_new_symbol(SYM_Temporary);
    
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), attr_list, i);
    free(attr_list);
    return result;
}

void attr_init(void) {
    symtab_add_builtin("Attributes", builtin_attributes);
    symtab_add_builtin("SetAttributes", builtin_set_attributes);
    symtab_add_builtin("ClearAttributes", builtin_clear_attributes);
    
    /* Seed the builtin attributes into the symbol table, once. This is a
     * starting value, not a floor: the bits live only in def->attributes, so
     * a later Unprotect / ClearAttributes is not silently undone (issue #103). */
    for (int i = 0; builtin_attrs[i].name != NULL; i++) {
        SymbolDef* def = symtab_get_def(builtin_attrs[i].name);
        if (def) def->attributes |= builtin_attrs[i].attrs;
    }
}
