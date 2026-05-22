#include "attr.h"
#include "symtab.h"
#include "eval.h"
#include "sym_names.h"
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
    {"Unevaluated", ATTR_HOLDALLCOMPLETE | ATTR_PROTECTED},
    {"Set", ATTR_HOLDFIRST | ATTR_PROTECTED},
    {"SetDelayed", ATTR_HOLDALL | ATTR_PROTECTED},
    {"Clear", ATTR_HOLDALL | ATTR_PROTECTED},
    {"AppendTo", ATTR_HOLDFIRST | ATTR_PROTECTED},
    {"PrependTo", ATTR_HOLDFIRST | ATTR_PROTECTED},
    {"Plus", ATTR_FLAT | ATTR_ORDERLESS | ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_ONEIDENTITY | ATTR_PROTECTED},
    {"Times", ATTR_FLAT | ATTR_ORDERLESS | ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_ONEIDENTITY | ATTR_PROTECTED},
    {"Divide", ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED},
    {"Subtract", ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED},
    {"Power", ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_ONEIDENTITY | ATTR_PROTECTED},
    {"Sqrt", ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED},
    {"Rational", ATTR_PROTECTED},
    {"Attributes", ATTR_HOLDALL | ATTR_PROTECTED},
    {"SetAttributes", ATTR_HOLDFIRST | ATTR_PROTECTED},
    {"ClearAttributes", ATTR_HOLDFIRST | ATTR_PROTECTED},
    {"OwnValues", ATTR_HOLDALL | ATTR_PROTECTED},
    {"DownValues", ATTR_HOLDALL | ATTR_PROTECTED},
    {"Out", ATTR_PROTECTED},
    {"Overflow", ATTR_PROTECTED},
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
    {"Condition", ATTR_HOLDREST | ATTR_PROTECTED},
    {"Rule", ATTR_PROTECTED},
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
    {"Boole", ATTR_LISTABLE | ATTR_PROTECTED},
    {"CompoundExpression", ATTR_HOLDALL | ATTR_PROTECTED},
    {"Table", ATTR_HOLDALL | ATTR_PROTECTED},
    {"Module", ATTR_HOLDALL | ATTR_PROTECTED},
    {"Block", ATTR_HOLDALL | ATTR_PROTECTED},
    {"With", ATTR_HOLDALL | ATTR_PROTECTED},
    {"Range", ATTR_LISTABLE | ATTR_PROTECTED},
    {"Return", ATTR_PROTECTED},
    {"Array", ATTR_PROTECTED},
    {"Take", ATTR_NHOLDREST | ATTR_PROTECTED},
    {"Drop", ATTR_NHOLDREST | ATTR_PROTECTED},
    {"Extract", ATTR_NHOLDREST | ATTR_PROTECTED},
    {"Flatten", ATTR_PROTECTED},
    {"Partition", ATTR_PROTECTED},
    {"RotateLeft", ATTR_PROTECTED},
    {"RotateRight", ATTR_PROTECTED},
    {"Reverse", ATTR_PROTECTED},
    {"Transpose", ATTR_PROTECTED},
    {"Mod", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Quotient", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"QuotientRemainder", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"GCD", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE | ATTR_FLAT | ATTR_ORDERLESS | ATTR_ONEIDENTITY},
    {"LCM", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE | ATTR_FLAT | ATTR_ORDERLESS | ATTR_ONEIDENTITY},
    {"PrimeQ", ATTR_PROTECTED | ATTR_LISTABLE},
    {"PrimePi", ATTR_PROTECTED | ATTR_LISTABLE},
    {"FactorInteger", ATTR_PROTECTED | ATTR_LISTABLE},
    {"NextPrime", ATTR_LISTABLE | ATTR_PROTECTED | ATTR_READPROTECTED},
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
    {"Exp", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"E", ATTR_PROTECTED},
    {"Floor", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Ceiling", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Round", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"IntegerPart", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"IntegerDigits", ATTR_PROTECTED | ATTR_LISTABLE},
    {"FractionalPart", ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE},
    {"Timing", ATTR_HOLDALL | ATTR_PROTECTED | ATTR_SEQUENCEHOLD},
    {"RepeatedTiming", ATTR_HOLDFIRST | ATTR_PROTECTED | ATTR_SEQUENCEHOLD},
    {"TimeConstrained", ATTR_HOLDALL | ATTR_PROTECTED},
    {"Dot", ATTR_FLAT | ATTR_ONEIDENTITY | ATTR_PROTECTED},
    {"Det", ATTR_PROTECTED},
    {"Cross", ATTR_PROTECTED},
    {"Norm", ATTR_PROTECTED},
    {NULL, 0}
};

uint32_t get_attributes(const char* symbol_name) {
    if (!symbol_name) return ATTR_NONE;
    
    uint32_t attrs = ATTR_NONE;

    // Check builtin_attrs table first for base attributes
    for (int i = 0; builtin_attrs[i].name != NULL; i++) {
        if (strcmp(builtin_attrs[i].name, symbol_name) == 0) {
            attrs = builtin_attrs[i].attrs;
            break;
        }
    }

    // Merge with dynamic attributes from symbol table
    SymbolDef* def = symtab_get_def(symbol_name);
    if (def) {
        attrs |= def->attributes;
    }

    return attrs;
}

void set_attributes(const char* symbol_name, uint32_t attrs) {
    if (!symbol_name) return;
    SymbolDef* def = symtab_get_def(symbol_name);
    if (def) {
        // If the symbol is Locked, we cannot change its attributes
        if (def->attributes & ATTR_LOCKED) return;
        if (def->attributes != attrs) {
            def->attributes = attrs;
            /* Attribute changes (Hold*, Listable, Flat, Orderless,
             * NumericFunction, Protected, OneIdentity, ...) all change
             * how the evaluator handles a head, so a cached evaluation
             * could be stale. Bump the global eval clock to invalidate. */
            eval_clock_bump();
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
    if (strcmp(name, "NHoldRest") == 0) return ATTR_NHOLDREST;
    if (strcmp(name, "Locked") == 0) return ATTR_LOCKED;
    if (strcmp(name, "ReadProtected") == 0) return ATTR_READPROTECTED;
    if (strcmp(name, "Temporary") == 0) return ATTR_TEMPORARY;
    if (strcmp(name, "SequenceHold") == 0) return ATTR_SEQUENCEHOLD;
    return ATTR_NONE;
}

static void remove_single_attribute(SymbolDef* def, Expr* attr_expr) {
    const char* attr_name = NULL;
    if (attr_expr->type == EXPR_SYMBOL) attr_name = attr_expr->data.symbol;
    else if (attr_expr->type == EXPR_STRING) attr_name = attr_expr->data.string;

    if (attr_name) {
        uint32_t bit = string_to_attribute(attr_name);
        if (bit != ATTR_NONE && (def->attributes & bit)) {
            def->attributes &= ~bit;
            /* Real attribute change -- invalidate cached evaluations. */
            eval_clock_bump();
        }
    }
}

static void add_single_attribute(SymbolDef* def, Expr* attr_expr) {
    const char* attr_name = NULL;
    if (attr_expr->type == EXPR_SYMBOL) attr_name = attr_expr->data.symbol;
    else if (attr_expr->type == EXPR_STRING) attr_name = attr_expr->data.string;

    if (attr_name) {
        uint32_t bit = string_to_attribute(attr_name);
        if (bit != ATTR_NONE && !(def->attributes & bit)) {
            def->attributes |= bit;
            /* Real attribute change -- invalidate cached evaluations. */
            eval_clock_bump();
        }
    }
}

static void set_attributes_for_symbol(Expr* sym_expr, Expr* attr_spec) {
    const char* sym_name = NULL;
    if (sym_expr->type == EXPR_SYMBOL) sym_name = sym_expr->data.symbol;
    else if (sym_expr->type == EXPR_STRING) sym_name = sym_expr->data.string;
    
    if (!sym_name) return;
    
    SymbolDef* def = symtab_get_def(sym_name);
    if (!def || (def->attributes & ATTR_LOCKED)) return;

    // Attribute spec can be a single attribute or a list of attributes
    if (attr_spec->type == EXPR_SYMBOL || attr_spec->type == EXPR_STRING) {
        add_single_attribute(def, attr_spec);
    } else if (attr_spec->type == EXPR_FUNCTION && 
               attr_spec->data.function.head->type == EXPR_SYMBOL &&
               attr_spec->data.function.head->data.symbol == SYM_List) {
        for (size_t i = 0; i < attr_spec->data.function.arg_count; i++) {
            add_single_attribute(def, attr_spec->data.function.args[i]);
        }
    }
}

static void clear_attributes_for_symbol(Expr* sym_expr, Expr* attr_spec) {
    const char* sym_name = NULL;
    if (sym_expr->type == EXPR_SYMBOL) sym_name = sym_expr->data.symbol;
    else if (sym_expr->type == EXPR_STRING) sym_name = sym_expr->data.string;

    if (!sym_name) return;

    SymbolDef* def = symtab_get_def(sym_name);
    if (!def || (def->attributes & ATTR_LOCKED)) return;

    // Attribute spec can be a single attribute or a list of attributes
    if (attr_spec->type == EXPR_SYMBOL || attr_spec->type == EXPR_STRING) {
        remove_single_attribute(def, attr_spec);
    } else if (attr_spec->type == EXPR_FUNCTION &&
               attr_spec->data.function.head->type == EXPR_SYMBOL &&
               attr_spec->data.function.head->data.symbol == SYM_List) {
        for (size_t i = 0; i < attr_spec->data.function.arg_count; i++) {
            remove_single_attribute(def, attr_spec->data.function.args[i]);
        }
    }
}

Expr* builtin_clear_attributes(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;

    Expr* sym_spec = res->data.function.args[0];
    Expr* attr_spec = res->data.function.args[1];

    // Symbol spec can be a single symbol/string or a list of them
    if (sym_spec->type == EXPR_SYMBOL || sym_spec->type == EXPR_STRING) {
        clear_attributes_for_symbol(sym_spec, attr_spec);
    } else if (sym_spec->type == EXPR_FUNCTION &&
               sym_spec->data.function.head->type == EXPR_SYMBOL &&
               sym_spec->data.function.head->data.symbol == SYM_List) {
        for (size_t i = 0; i < sym_spec->data.function.arg_count; i++) {
            clear_attributes_for_symbol(sym_spec->data.function.args[i], attr_spec);
        }
    }

    return expr_new_symbol("Null");
}

Expr* builtin_set_attributes(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    
    Expr* sym_spec = res->data.function.args[0];
    Expr* attr_spec = res->data.function.args[1];
    
    // Symbol spec can be a single symbol/string or a list of them
    if (sym_spec->type == EXPR_SYMBOL || sym_spec->type == EXPR_STRING) {
        set_attributes_for_symbol(sym_spec, attr_spec);
    } else if (sym_spec->type == EXPR_FUNCTION && 
               sym_spec->data.function.head->type == EXPR_SYMBOL &&
               sym_spec->data.function.head->data.symbol == SYM_List) {
        for (size_t i = 0; i < sym_spec->data.function.arg_count; i++) {
            set_attributes_for_symbol(sym_spec->data.function.args[i], attr_spec);
        }
    }
    
    return expr_new_symbol("Null");
}

Expr* builtin_attributes(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    const char* sym_name = NULL;
    if (arg->type == EXPR_SYMBOL) sym_name = arg->data.symbol;
    else if (arg->type == EXPR_STRING) sym_name = arg->data.string;
    
    if (!sym_name) return NULL;
    
    uint32_t attrs = get_attributes(sym_name);
    
    // Count attributes
    size_t count = 0;
    if (attrs & ATTR_HOLDFIRST) count++;
    if (attrs & ATTR_HOLDREST) count++;
    if (attrs & ATTR_HOLDALLCOMPLETE) count++;
    if (attrs & ATTR_FLAT) count++;
    if (attrs & ATTR_ORDERLESS) count++;
    if (attrs & ATTR_LISTABLE) count++;
    if (attrs & ATTR_NUMERICFUNCTION) count++;
    if (attrs & ATTR_PROTECTED) count++;
    if (attrs & ATTR_ONEIDENTITY) count++;
    if (attrs & ATTR_NHOLDREST) count++;
    if (attrs & ATTR_LOCKED) count++;
    if (attrs & ATTR_READPROTECTED) count++;
    if (attrs & ATTR_TEMPORARY) count++;
    if (attrs & ATTR_SEQUENCEHOLD) count++;

    Expr** attr_list = malloc(sizeof(Expr*) * count);
    size_t i = 0;
    if (attrs & ATTR_FLAT) attr_list[i++] = expr_new_symbol("Flat");
    if ((attrs & ATTR_HOLDALL) == ATTR_HOLDALL) {
        attr_list[i++] = expr_new_symbol("HoldAll");
    } else {
        if (attrs & ATTR_HOLDFIRST) attr_list[i++] = expr_new_symbol("HoldFirst");
        if (attrs & ATTR_HOLDREST) attr_list[i++] = expr_new_symbol("HoldRest");
    }
    if (attrs & ATTR_HOLDALLCOMPLETE) attr_list[i++] = expr_new_symbol("HoldAllComplete");
    if (attrs & ATTR_LISTABLE) attr_list[i++] = expr_new_symbol("Listable");
    if (attrs & ATTR_LOCKED) attr_list[i++] = expr_new_symbol("Locked");
    if (attrs & ATTR_NUMERICFUNCTION) attr_list[i++] = expr_new_symbol("NumericFunction");
    if (attrs & ATTR_ONEIDENTITY) attr_list[i++] = expr_new_symbol("OneIdentity");
    if (attrs & ATTR_NHOLDREST) attr_list[i++] = expr_new_symbol("NHoldRest");
    if (attrs & ATTR_ORDERLESS) attr_list[i++] = expr_new_symbol("Orderless");
    if (attrs & ATTR_PROTECTED) attr_list[i++] = expr_new_symbol("Protected");
    if (attrs & ATTR_READPROTECTED) attr_list[i++] = expr_new_symbol("ReadProtected");
    if (attrs & ATTR_SEQUENCEHOLD) attr_list[i++] = expr_new_symbol("SequenceHold");
    if (attrs & ATTR_TEMPORARY) attr_list[i++] = expr_new_symbol("Temporary");
    
    Expr* result = expr_new_function(expr_new_symbol("List"), attr_list, i);
    free(attr_list);
    return result;
}

void attr_init(void) {
    symtab_add_builtin("Attributes", builtin_attributes);
    symtab_add_builtin("SetAttributes", builtin_set_attributes);
    symtab_add_builtin("ClearAttributes", builtin_clear_attributes);
    
    // Initialize builtin attributes in symtab
    for (int i = 0; builtin_attrs[i].name != NULL; i++) {
        SymbolDef* def = symtab_get_def(builtin_attrs[i].name);
        if (def) {
            def->attributes = builtin_attrs[i].attrs;
        }
    }
}
