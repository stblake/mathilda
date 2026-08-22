/* Mathilda — Compile[] engine: emitter-internal shared surface.
 *
 * The emitter (compile.c) and its sibling emit modules — compile_mgd.c,
 * compile_assoc.c, compile_infer.c and the compile_emit_*.c category files —
 * all build ONE program through the compile-time `Ctx` builder.  `Ctx` / `Val`,
 * the emitter constants, and the register-allocation / instruction-emit /
 * lowering primitives live here so those files can share them.
 *
 * NONE of this is visible to the VM (compile_vm.c) or to the optimiser /
 * disassembler, which include only compile_internal.h — the finished-program
 * representation.  Keeping the two headers apart is what stops the runtime ever
 * depending on the emitter's internal builder types.
 *
 * Not part of the public API — clients hold an opaque CompiledProgram and go
 * through compile.h. */
#ifndef MATHILDA_COMPILE_EMIT_H
#define MATHILDA_COMPILE_EMIT_H

#include "compile_internal.h"   /* Slot / Instr / CompiledProgram / PartSpec / AssocSpec / MgdConst */
#include "../ndarray.h"         /* NDUnaryKernel — in nd_unary_elem's signature */

struct CompiledFunction;        /* compiled_callee's return; full type in compiled_function.h */

/* ------------------------------------------------------------------ *
 *  Emitter constants                                                  *
 * ------------------------------------------------------------------ */
/* A compiled callee bigger than this is CALLed rather than pasted in: past this
 * size inlining costs more in code size and compile time than the call costs at
 * run time. */
#define INLINE_MAX_INSTRS 32

#define CSE_MAX       16     /* reserved CSE registers; bounds frame growth */
#define CSE_OCC_MAX  256     /* occurrence map entries */

/* Depth of the lexical scope stack.  One entry per live binding, and the
 * bindings now come from six places at once — Sum/Product iterators, Do
 * iterators, With/Module locals, inlined CompiledFunction parameters, fusion
 * leaves, and (since the functional heads landed) every lambda parameter of
 * every Nest/Fold/Map in the body.  Nesting three of those exhausted 16, and a
 * scope overflow bails the WHOLE body, so the array is sized for nesting rather
 * than for a single construct.  It lives on the stack-allocated Ctx, so the
 * cost is bytes at compile time and nothing at run time. */
#define CTX_MAX_SCOPE 32

/* Live Slot[] bindings — the parameters of a `Function[body]` being inlined.
 * Also bounds a lambda's parameter count and hence emit_apply's argument
 * count, so the three stay consistent by construction. */
#define FN_MAX_PARAMS 8

/* Rungs in a first-match-wins conditional ladder (Which / Switch / Piecewise).
 * Far past any real ladder; a longer one bails to the interpreter, which is
 * always correct, just slower. */
#define COND_CHAIN_MAX 256

/* ------------------------------------------------------------------ *
 *  Argument name -> index map (interned-pointer hash)                 *
 * ------------------------------------------------------------------ */
typedef struct { const char** key; int* val; size_t cap; } NameMap;

/* A lowered value: the register that holds it, whether that register is a temp
 * the caller must release, its compile type, and whether an array value was
 * CONSTRUCTED by the body (so the interpreter would answer a List). */
typedef struct { int reg; bool tmp; CompileType type; bool built; } Val;

/* ------------------------------------------------------------------ *
 *  Compiler context (the program builder)                             *
 * ------------------------------------------------------------------ */
typedef struct {
    Instr* code; size_t n, cap;
    int nlocals;        /* registers [0,nlocals) are args/locals */
    int temp_top;       /* current scalar temp-stack height      */
    int maxreg;         /* high-water scalar register count      */
    int arr_top;        /* current array temp-stack height       */
    int arr_max;        /* high-water array register count       */
    const CompileType* arg_types;
    NameMap map;
    unsigned char* argdep;
    /* lexically-scoped loop variables (Sum/Product/Nest), innermost last */
    struct { const char* name; int reg; CompileType type; bool built; } scope[CTX_MAX_SCOPE];
    int nscope;

    /* Live Slot[] bindings for an inlined `Function[body]` (`#`, `#1`, `#2`).
     *
     * Slot[k] is an EXPR_FUNCTION, not a symbol, so it cannot ride on scope[],
     * which matches interned symbol POINTERS.  One flat frame is not an
     * approximation: substitute_slots (src/purefunc.c:76) deliberately does not
     * recurse into a nested Function, so only the innermost frame is ever
     * visible, and emit_apply saves/restores nslot around each body. */
    struct { int reg; CompileType type; } slot[FN_MAX_PARAMS];
    int nslot;
    bool ok;
    unsigned flags;     /* COMPILE_FOLD_GLOBALS, ... */
    int inlining;       /* >0 while lowering an inlined CompiledFunction body */
    int fusing;         /* >0 while lowering the body of a fused elementwise loop */
    int tile_top, tile_max;   /* strip-mining tile registers (TILE_VREG bank) */
    bool vector_mode;   /* inside a strip-mined loop: tile-valued ops are legal */
    bool array_args;    /* any declared argument is an array: gates fusion probing */

    /* General-Part subscript lists (see PartSpec).  Held here while emitting and
     * handed to the program at finalize; freed here if the compile bails. */
    PartSpec** parts;
    int        nparts, parts_cap;

    /* Association read-op specs (see AssocSpec).  Same lifetime as `parts`:
     * built while emitting, handed to the program at finalize, freed here on a
     * bail. */
    AssocSpec** assocs;
    int         nassocs, assocs_cap;

    /* B4 higher-order transform callees (see AssocCalleeSpec), same lifetime. */
    AssocCalleeSpec** callees;
    int         ncallees, callees_cap;

    /* Expr-level CSE (see cse_plan).  A chosen repeated subtree is computed once
     * into a register reserved below the temp stack, so no temp allocation and
     * none of the explicit temp_top resets in the loop lowerings can reach it. */
    struct { const Expr* node; int idx; } cse_occ[CSE_OCC_MAX];
    size_t      cse_nocc;
    int         cse_reg[CSE_MAX];
    CompileType cse_type[CSE_MAX];
    bool        cse_ready[CSE_MAX];
    int         ncse;

    /* Diagnostics: the INNERMOST subexpression `emit` could not lower.  A bail is
     * otherwise invisible — the caller just interprets, an order of magnitude
     * slower, and a body one head outside the subset costs the whole body.  Set
     * in the emit wrapper (first writer wins, so recursion unwinding leaves the
     * deepest cause) and rolled back with speculative lowering.  Borrowed: the
     * body outlives the compile. */
    const Expr* bail_node;

    /* Arbitrary precision (see emit_mgd).  Both zero/false for a machine program,
     * in which case none of the managed machinery below is ever touched and the
     * machine emit path (emit_node) is used unchanged. */
    long prec_bits;         /* MPFR working precision in bits, or 0 (machine)   */
    bool allow_bigint;      /* "BigIntegers" -> True                            */
    int  mgd_top, mgd_max;  /* managed-bank alloc (bump-only; never popped, so   *
                             * a physical managed register has ONE container type)*/
    CompileType* mgd_type;  /* [mgd_max] container type per managed vreg index    */
    int  mgd_type_cap;
    MgdConst* mgd_consts;   /* program-owned literal containers; handed to program*/
    int  nmgd_consts, mgd_consts_cap;
} Ctx;

/* ------------------------------------------------------------------ *
 *  Compile-time function values, loop specs, integer-closed heads     *
 * ------------------------------------------------------------------ *
 * Shared between the core emitter, the inference pass (compile_infer.c) and the
 * category emitters: FnSpec describes a resolved `Function[...]` / head value;
 * LoopSpec a counted iteration; IntClosed an integer-closed head. */
typedef enum {
    FN_LAMBDA,    /* Function[u, body] / Function[{u, ...}, body] */
    FN_SLOTS,     /* Function[body] with #, #1, #2, ...           */
    FN_HEAD,      /* a bare head: Sin, Plus, or a CompiledFunction-valued symbol */
    FN_IDENTITY,  /* Identity                                     */
    FN_COMPOSE    /* Composition[f1, ..., fn]                     */
} FnKind;

typedef struct {
    FnKind      kind;
    int         nparams;                /* -1 == accepts any arity */
    const char* pname[FN_MAX_PARAMS];   /* FN_LAMBDA: interned parameter names */
    const Expr* body;                   /* FN_LAMBDA / FN_SLOTS */
    const char* head;                   /* FN_HEAD: interned head name */
    /* The `f` node as written.  FN_COMPOSE walks it; every kind uses it to
     * blame a bail on the user's own tree rather than on scaffolding. */
    const Expr* fexpr;
} FnSpec;

typedef struct {
    const char* var;    /* iterator symbol, or NULL for the count-only forms */
    const Expr* lo;     /* NULL => literal 1 */
    const Expr* hi;
    long long   di;     /* step; nonzero */
} LoopSpec;

/* One rung of a first-match-wins conditional ladder (Which / Switch / Piecewise):
 * a guard and the value it yields.  `key == NULL` marks the unconditional default
 * rung (Piecewise's default, Which's `True`, Switch's `_`), which is always last.
 * For Switch, `key` is a literal form and the guard is `discriminant == key`; for
 * Which/Piecewise `key` is the boolean condition itself. */
typedef struct { const Expr* key; const Expr* value; } Rung;

typedef struct { const char* name; size_t arity; uint16_t op; } IntClosed;

/* Snapshot of the builder state for speculative lowering: emit_mark captures it,
 * emit_rollback restores it so a failed probe (fusion, a compiled callee) leaves
 * no half-emitted instructions or leaked temps behind. */
typedef struct {
    size_t n; int temp_top, maxreg, arr_top, arr_max, nscope; bool ok;
    int tile_top, tile_max; bool vector_mode;
    const Expr* bail_node;
} EmitMark;

/* Array-valued leaves of an elementwise chain being probed for strip-mined
 * fusion (fuse_collect gathers them; fuse_listable / try_fuse drive the probe). */
#define FUSE_MAX_LEAVES 8
typedef struct {
    const char* name[FUSE_MAX_LEAVES];   /* interned symbol pointer */
    int         reg[FUSE_MAX_LEAVES];    /* the array-valued register */
    CompileType type[FUSE_MAX_LEAVES];   /* the ARRAY type (element + rank) */
    bool        built[FUSE_MAX_LEAVES];  /* leaf constructed here, not an argument */
    int         n;
} FuseLeaves;

/* ------------------------------------------------------------------ *
 *  Promoted emit-core primitives (defined in compile.c)               *
 * ------------------------------------------------------------------ *
 * Small building blocks the extracted emit modules share with the core.  Each
 * is `static` no longer only because a sibling file calls it; none is part of
 * the public API. */
int  arg_find(const Ctx* c, const char* nm);
int  scope_find(const Ctx* c, const char* nm, CompileType* type, bool* built);
bool literal(const Expr* e, Slot* imm, CompileType* type);
CompileType num_common(CompileType a, CompileType b);
void ins(Ctx* c, uint16_t op, uint32_t dst, uint32_t a, uint32_t b, Slot imm);
void ins_f(Ctx* c, uint16_t op, uint16_t flags, uint32_t dst, uint32_t a, uint32_t b, Slot imm);
int  alloc_temp(Ctx* c);
int  alloc_arr(Ctx* c);
void pop_tmp(Ctx* c, Val v);
void free_if_tmp(Ctx* c, Val v);
void coerce(Ctx* c, Val* v, CompileType target);
Val  emit_const(Ctx* c, Slot imm, CompileType type);
bool emit(Ctx* c, const Expr* e, Val* out);
bool infer_type(Ctx* c, const Expr* e, CompileType* out);
struct CompiledProgram* compile_value_callee(Ctx* c, const Expr* fn, CompileType in);
bool global_const(const Ctx* c, const char* nm, Slot* imm, CompileType* type);
const struct CompiledFunction* compiled_callee(const Ctx* c, const char* nm);
bool ct_is_elem(CompileType t);
bool named_const(const char* nm, double* out);
bool unary_math(const char* h, uint16_t* op_r, uint16_t* op_c);
CompileType nd_unary_elem(const NDUnaryKernel* k, CompileType ea);
bool fn_resolve(const Expr* f, int want_arity, FnSpec* out);
int  fn_slot_index(const Ctx* c, Expr* const* A, size_t na);
Expr* fn_head_call(const char* head, int n);
const char* fn_placeholder(int i);
bool loop_spec_parse(const Expr* spec, LoopSpec* out);
bool loop_spec_int_bounds(Ctx* c, const LoopSpec* s);
/* Walk a Which / Switch / Piecewise call into the live rungs of its first-match
 * ladder and the ladder's result type, dropping the clauses the interpreter drops
 * (a literal-False guard) and stopping at the one it makes final (a literal-True
 * guard, Switch's `_`).  Shared by emit_ctrl (which needs the rungs) and infer_type
 * (which passes rungs=NULL and reads only *rt) so the declared type and the emitted
 * code cannot drift.  Returns false — with c->ok cleared — for anything outside the
 * compilable subset. */
bool collect_ladder(Ctx* c, const char* h, Expr** A, size_t na,
                    Rung* rungs, int* nr_out, bool* has_default_out, CompileType* rt_out);
bool range_spec(Ctx* c, Expr* const* A, size_t na, LoopSpec* out);
bool subscript_is_literal_spec(const Expr* e);
bool part_is_scalar_indexed(Ctx* c, CompileType at, const Expr* const* A, size_t na);
bool const_array_shape(Ctx* c, const Expr* const* A, int* rank_out, CompileType* elem_out);
CompileType nd_binary_elem(const NDBinaryKernel* k, CompileType ea, CompileType es);

/* Register/value machinery + scalar/vector/array op emitters + guards, shared
 * with the per-category emitters (compile_emit_*.c). */
int  alloc_tile(Ctx* c);
bool val_is_tile(Val v);
bool reg_is_owned_arr(int r);
bool scalar_only(Ctx* c, CompileType a, CompileType b, CompileType r);
bool expr_subtree_of(const Expr* root, const Expr* node);
bool ctx_own_partspec(Ctx* c, PartSpec* p);
int  cse_lookup(const Ctx* c, const Expr* e);
Val  binop(Ctx* c, uint16_t op, Val a, Val b, CompileType rtype);
Val  unop(Ctx* c, uint16_t op, Val a, CompileType rtype);
Val  vsplat(Ctx* c, Val v);
Val  vec_unop(Ctx* c, uint16_t op, Val a, CompileType rtype, Slot imm);
Val  kern_unop(Ctx* c, uint16_t op, Val a, CompileType rt, const void* fn);
Val  kern_binop(Ctx* c, uint16_t op, Val a, Val b, CompileType rt, const void* fn);
Val  arr_op(Ctx* c, uint16_t op, Val a, Val b, CompileType rt, Slot imm);
Val  arr_noop_val(void);
void arr_prep(Ctx* c, Val* v, CompileType elem);
Val  arr_real_const(Ctx* c, double x);
Val  arr_ew(Ctx* c, Val a, Val b, CompileType rt, bool is_plus);
bool emit_flat_index(Ctx* c, Val arr, const Expr* const* A, size_t na, int* idx_out);
PartSpec* emit_partspec(Ctx* c, const Expr* const* A, size_t na, int* base_out);
bool emit_arr_unary(Ctx* c, const char* head, Val a, Val* out);
bool emit_unary_math(Ctx* c, const char* head, const Expr* arg, uint16_t op_r, uint16_t op_c, Val* out);
bool emit_loop_bound(Ctx* c, const Expr* b, int dst);
bool try_kernel(Ctx* c, const char* h, Expr** A, size_t na, Val* out);
EmitMark emit_mark(const Ctx* c);
void emit_rollback(Ctx* c, EmitMark m);
bool emit_apply(Ctx* c, const FnSpec* s, const Val* argv, int n, Val* out);
void emit_max_guard(Ctx* c, int reg, long long lim);
bool emit_iter_count(Ctx* c, const LoopSpec* s, int rlo, int dst);
bool stmt_valued_head(const Expr* e);
uint16_t emit_sameq_op(CompileType t);
void emit_nonneg_guard(Ctx* c, int reg);
void emit_nonzero_guard(Ctx* c, int reg);
uint16_t a_load_op(CompileType elem);
uint16_t a_store_op(CompileType elem);
bool any_array_arg(Ctx* c, Expr** A, size_t na);
bool fuse_collect(Ctx* c, const Expr* e, FuseLeaves* L);
int  fuse_listable(Ctx* c, const Expr* e);

/* The seven per-category emitters (compile_emit_*.c).  Each returns 1 if `h` is
 * one of ITS heads (then *out is valid iff c->ok), 0 if not handled — emit_node
 * calls them in sequence and falls through to the general kernel dispatch. */
int emit_arith (Ctx* c, const char* h, const Expr* e, Expr** A, size_t na, Val* out);
int emit_mathfn(Ctx* c, const char* h, const Expr* e, Expr** A, size_t na, Val* out);
int emit_int   (Ctx* c, const char* h, const Expr* e, Expr** A, size_t na, Val* out);
int emit_array (Ctx* c, const char* h, const Expr* e, Expr** A, size_t na, Val* out);
int emit_logic (Ctx* c, const char* h, const Expr* e, Expr** A, size_t na, Val* out);
int emit_ctrl  (Ctx* c, const char* h, const Expr* e, Expr** A, size_t na, Val* out);
int emit_iter  (Ctx* c, const char* h, const Expr* e, Expr** A, size_t na, Val* out);

/* ------------------------------------------------------------------ *
 *  Module entry points (each lowered in its own .c)                   *
 * ------------------------------------------------------------------ */
bool emit_mgd(Ctx* c, const Expr* e, Val* out);   /* compile_mgd.c — managed scalars */
/* Association B1–B5 (compile_assoc.c).  try_emit_assoc handles an association
 * read op or returns false so the ordinary lowering takes the node;
 * try_infer_assoc is its inference-pass twin. */
bool try_emit_assoc(Ctx* c, const char* h, Expr** A, size_t na, Val* out);
bool try_infer_assoc(Ctx* c, const char* h, Expr** A, size_t na, CompileType* out);
/* Type-inference engine (compile_infer.c): infer_type is the entry; the rest
 * type the functional / integer-closed heads for emit_node's result register. */
bool infer_apply(Ctx* c, const FnSpec* s, const CompileType* argt, int n, CompileType* out);
int  nest_fixed_type(Ctx* c, const FnSpec* s, CompileType t0);
int  accum_fixed_type(Ctx* c, const FnSpec* s, CompileType t0, const CompileType* rest, int nrest);
CompileType vec_elem_type(Ctx* c, const Expr* e);
bool fp_opts(Expr* const* A, size_t na, size_t start, const Expr** max_out, const Expr** same_out);
const IntClosed* int_closed_head(const char* h, size_t na);
bool int_only_head(const char* h, size_t na, bool* pred);
bool all_args_int(Ctx* c, Expr* const* A, size_t na);

#endif /* MATHILDA_COMPILE_EMIT_H */
