/*
 * opform.c -- generic operator (curried) forms h[o1..ok][x].  See opform.h.
 *
 * The table is a small array searched by interned-name pointer.  It is only
 * consulted for a composite head h[...] applied to exactly one argument that
 * no earlier rule in the evaluator's composite-head chain claimed, so a
 * linear scan over a few dozen pointers is far below the cost of the
 * evaluation it enables.
 */
#include "opform.h"
#include "eval.h"
#include "sym_intern.h"
#include "symtab.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const char* head;      /* interned symbol name */
    size_t      min_ops;   /* fewest operator arguments accepted */
    size_t      max_ops;   /* most operator arguments accepted */
    size_t      insert_at; /* 0-based slot of the data argument x */
} OpFormSpec;

#define OPFORM_MAX 256
static OpFormSpec g_opforms[OPFORM_MAX];
static size_t     g_opform_count = 0;

static const OpFormSpec* opform_find(const char* interned) {
    for (size_t i = 0; i < g_opform_count; i++)
        if (g_opforms[i].head == interned) return &g_opforms[i];
    return NULL;
}

/* Append "h[o1][x] is h[x, o1]" usage lines to the head's docstring (unless it
 * already describes an operator form), so ?h always documents the curried
 * call it accepts. */
static void opform_document(const OpFormSpec* row) {
    const char* old = symtab_get_docstring(row->head);
    if (old && (strstr(old, "operator form") || strstr(old, "Operator form"))) return;
    char buf[1024];
    size_t len = 0;
    len += (size_t)snprintf(buf + len, sizeof buf - len, "%sOperator form:", old ? "\n" : "");
    for (size_t k = row->min_ops; k <= row->max_ops && len < sizeof buf - 128; k++) {
        char ops[64] = "", full[96] = "";
        size_t ol = 0, fl = 0;
        for (size_t i = 0; i < k; i++)
            ol += (size_t)snprintf(ops + ol, sizeof ops - ol, "%sa%zu", i ? ", " : "", i + 1);
        for (size_t i = 0, j = 0; i <= k; i++) {
            if (i == row->insert_at) fl += (size_t)snprintf(full + fl, sizeof full - fl, "%sx", i ? ", " : "");
            else { j++; fl += (size_t)snprintf(full + fl, sizeof full - fl, "%sa%zu", i ? ", " : "", j); }
        }
        len += (size_t)snprintf(buf + len, sizeof buf - len, "%s %s[%s][x] is %s[%s]",
                                k > row->min_ops ? ";" : "", row->head, ops, row->head, full);
    }
    if (len >= sizeof buf) len = sizeof buf - 1;
    buf[len] = '\0';
    size_t ol = old ? strlen(old) : 0;
    char* doc = malloc(ol + len + 2);
    if (!doc) return;
    memcpy(doc, old ? old : "", ol);
    memcpy(doc + ol, buf, len);
    doc[ol + len] = '.';
    doc[ol + len + 1] = '\0';
    symtab_set_docstring(row->head, doc);
    free(doc);
}

void opform_register(const char* head, size_t min_ops, size_t max_ops, size_t insert_at) {
    if (!head || min_ops > max_ops || insert_at > min_ops) return;
    const char* h = intern_symbol(head);
    OpFormSpec* row = (OpFormSpec*)opform_find(h);
    if (!row) {
        if (g_opform_count >= OPFORM_MAX) return;
        row = &g_opforms[g_opform_count++];
    }
    row->head = h;
    row->min_ops = min_ops;
    row->max_ops = max_ops;
    row->insert_at = insert_at;
    opform_document(row);
}

static const OpFormSpec* opform_spec_for(const Expr* head) {
    if (!head || head->type != EXPR_FUNCTION) return NULL;
    const Expr* hh = head->data.function.head;
    if (!hh || hh->type != EXPR_SYMBOL) return NULL;
    const OpFormSpec* spec = opform_find(hh->data.symbol.name);
    if (!spec) return NULL;
    size_t k = head->data.function.arg_count;
    if (k < spec->min_ops || k > spec->max_ops) return NULL;
    return spec;
}

bool opform_matches(const Expr* head) {
    return opform_spec_for(head) != NULL;
}

Expr* opform_apply(const Expr* head, const Expr* x) {
    const OpFormSpec* spec = opform_spec_for(head);
    if (!spec || !x) return NULL;

    size_t k = head->data.function.arg_count;
    size_t n = k + 1;
    Expr* stack_args[8];
    Expr** args = (n <= 8) ? stack_args : malloc(sizeof(Expr*) * n);
    for (size_t i = 0, j = 0; i < n; i++) {
        if (i == spec->insert_at) args[i] = expr_copy((Expr*)x);
        else                      args[i] = expr_copy(head->data.function.args[j++]);
    }
    Expr* call = expr_new_function(expr_copy(head->data.function.head), args, n);
    if (args != stack_args) free(args);

    Expr* result = evaluate(call);
    /* The head declined these arguments: keep the curried form, as
     * Mathematica does, instead of exposing the rewritten call. */
    if (result && expr_eq(result, call)) {
        expr_free(result);
        result = NULL;
    }
    expr_free(call);
    return result;
}

void opform_init(void) {
    /* x first: h[o][x] == h[x, o].  Every row below was checked against
     * Mathematica 15; heads without an operator form there (Take, Drop, ...)
     * are deliberately absent. */
    static const char* const DATA_FIRST_1[] = {
        "Select", "SelectFirst", "AllTrue", "AnyTrue", "NoneTrue",
        "KeyTake", "KeyDrop", "KeySelect", "KeySortBy",
        "KeyExistsQ", "KeyMemberQ", "KeyFreeQ", "FreeQ",
        "GroupBy", "CountsBy", "Merge", "DeleteDuplicatesBy",
        "ReplaceAll", "Replace", "ReplacePart",
        "Append", "Prepend", "Delete",
        "TakeLargest", "TakeSmallest",
    };
    for (size_t i = 0; i < sizeof(DATA_FIRST_1) / sizeof(DATA_FIRST_1[0]); i++)
        opform_register(DATA_FIRST_1[i], 1, 1, 0);

    opform_register("Lookup", 1, 1, 0);          /* Lookup[k][a]; WL has no Lookup[k, d][a] */
    opform_register("Insert", 2, 2, 0);          /* Insert[e, n][x]  == Insert[x, e, n] */
    opform_register("TakeLargestBy", 2, 2, 0);   /* TakeLargestBy[f, n][x] */
    opform_register("TakeSmallestBy", 2, 2, 0);

    /* x last: h[f][x] == h[f, x]. */
    static const char* const DATA_SECOND_1[] = {
        "Map", "Apply", "KeyMap", "KeyValueMap", "AssociationMap",
    };
    for (size_t i = 0; i < sizeof(DATA_SECOND_1) / sizeof(DATA_SECOND_1[0]); i++)
        opform_register(DATA_SECOND_1[i], 1, 1, 1);
}
