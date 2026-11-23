#ifndef MATHILDA_COMPILE_HANDLE_OPS_H
#define MATHILDA_COMPILE_HANDLE_OPS_H

/* ---------------------------------------------------------------------------
 * CompileHandleOps -- one row per Compile[] handle kind.
 *
 * The compiler's value types are machine scalars plus a few HANDLE kinds whose
 * register holds an Expr* in the `.arr` slot: packed arrays, read-only
 * Associations, LinkedLists and PVectors.  Each used to be hard-coded at every
 * site that cared (about 90): "is this a handle", "reject it as a scalar
 * operand", "its display name", "retain it instead of copying it", "bind it as
 * an argument", "hand it back as a result".  This table is the Compile half of
 * the project's ds_ops registry (src/ds_ops.h is the interpreter half): a kind
 * states those facts once, and the shared sites ask the row, so a new handle
 * kind is one row here rather than an edit at every site.
 *
 * What stays out of the table on purpose: the lowering of a kind's own
 * operations (Append, Lookup, Part, ...) is genuinely kind-specific and lives in
 * the kind's own module; the CompileType numbering and every opcode are
 * unchanged (a row only describes the encoding window the macros in compile.h
 * already define).
 *
 * A new kind is one more row in the table in handle_ops.c; lookup scans the
 * handful of windows (compile-time and call-boundary code only, never a VM
 * dispatch).
 * -------------------------------------------------------------------------- */

#include "compile.h"
#include <stdbool.h>
#include <stdint.h>

enum {
    /* Takes part in numeric broadcasting (num_common); the other kinds are never
     * a numeric operand. */
    CHF_NUMERIC_OPERAND = 1u << 0,
    /* VALUE semantics: an update yields a new value and the old one is kept, so
     * the kind may not be rebound by a With constant, joined by an If, or flow
     * through a function application (whose lowerings MOVE registers). */
    CHF_VALUE_SEMANTICS = 1u << 1,
    /* May be a scope (Module/Block/With) local and the scope's result, carried
     * in the owned array bank. */
    CHF_SCOPE_LOCAL     = 1u << 2
};

typedef struct CompileHandleOps {
    const char* name;        /* ct_name() display text; also CompileDiagnostics */
    CompileType lo;          /* first encoding of the window                    */
    int         span;        /* width of the window, in encodings (multiple of 4) */
    unsigned    flags;       /* CHF_*                                           */
    /* Opcode that takes one more reference to a handle register (OP_LL_RETAIN /
     * OP_PV_RETAIN), or 0 for a kind with no shared ownership (a borrowed
     * array / association is copied or refused instead). */
    uint16_t    retain_op;
    /* Element kinds the parser may leave open (declared without one: element
     * type CT_BOOL) and the compiler then searches: how many candidates (0 =
     * the kind never has an open element), and the bail reason when the open
     * element reaches the emitter. */
    int         open_choices;
    const char* open_bail;
    /* Marshalling boundary, kind-specific only where the checks differ.
     * borrow_arg: validate a boxed Expr as an argument of type `t` (BORROWED,
     * stored as-is by the caller); NULL = the kind has its own site (arrays pack
     * nested Lists there).  load_arg: bind a boxed argument of declared type
     * `want` to a register; sets *arr to the handle, false = decline. */
    bool (*borrow_arg)(const Expr* e, CompileType t);
    bool (*load_arg)(const CompileValue* v, CompileType want, Expr** arr);
} CompileHandleOps;

/* The row whose window contains `t`, or NULL (scalars, managed types, errors).
 * CT_IS_HANDLE(t) is `ct_is_handle(t)`. */
const CompileHandleOps* ct_handle_ops(CompileType t);

/* Element type of a handle encoding (the low two bits of its offset into the
 * window).  `t` must satisfy ct_is_handle. */
CompileType ct_handle_elem(CompileType t);

/* Display name of a handle type, or NULL when `t` is not a handle. */
const char* ct_handle_name(CompileType t);

/* A handle that is not a numeric operand (an Association, LinkedList or
 * PVector): what num_common refuses before it looks at array ranks. */
bool ct_is_opaque_handle(CompileType t);

/* Flag queries: false for a non-handle. */
bool ct_has_value_semantics(CompileType t);   /* CHF_VALUE_SEMANTICS */
bool ct_is_scope_local(CompileType t);        /* CHF_SCOPE_LOCAL     */

/* The retain opcode of a handle type (OP_LL_RETAIN / OP_PV_RETAIN), 0 when the
 * kind has none or `t` is not a handle. */
uint16_t ct_retain_op(CompileType t);

#endif /* MATHILDA_COMPILE_HANDLE_OPS_H */
