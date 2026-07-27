/* Mathilda — Compile[]: bytecode disassembler (the engine behind CompilePrint).
 *
 * Code generation is otherwise invisible.  A body either compiles or bails, and
 * the only observable facts are "it ran fast" and an instruction count from
 * CompileDiagnostics — which cannot answer the questions that actually come up:
 * which twelve instructions the optimiser left, whether the constants got folded
 * into them (K_BINK) or are still costing a CONST each, whether an array chain
 * fused into a strip-mined loop or is delegating per operation, whether OP_APAR
 * was emitted so the map will really fan out, and which machine kernel a special
 * function lowered to.
 *
 * Rendering rules that keep the output trustworthy:
 *   - Operand SHAPE always comes from `compile_op_kind[]`, never from a private
 *     switch on the opcode.  That is the same table the optimiser and patch_reg
 *     use, so a new opcode cannot be decoded one way here and another way there.
 *   - No raw addresses.  Kernel pointers are reverse-resolved to their symbol
 *     names through the symbol table, callee programs and parallel loops are
 *     numbered.  The output is therefore deterministic, diffable and testable.
 *   - An unrecognised opcode degrades to its raw fields rather than printing
 *     something confidently wrong.
 *
 * Registers are named by BANK — R scalar, V array handle, T strip-mining tile —
 * and numbered by frame slot, so a number means the same thing everywhere,
 * including in a branch target.
 */
#include "compile_internal.h"
#include "compiled_function.h"
#include "../symtab.h"
#include "../ndarray.h"
#include "../print.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <complex.h>

/* The fourth expansion of OPLIST (after the enum, the VM's jump table and the
 * optimiser's kind table): opcode -> name.  Generated from the same list, so a
 * new opcode cannot be added without one. */
#define X(name, kind) #name,
static const char* const op_name[OP__COUNT] = { OPLIST };
#undef X

/* ------------------------------------------------------------------ *
 *  Output buffer                                                      *
 * ------------------------------------------------------------------ */
/* Same growable-string shape as LBuf in print_latex.c, with one change: the
 * formatted append is sized with a two-pass vsnprintf instead of a fixed
 * scratch array, because a Part subscript or a printed body has no length
 * bound and a silently truncated disassembly is worse than none. */
typedef struct { char* s; size_t len, cap; } DBuf;

static void db_init(DBuf* b) {
    b->cap = 1024; b->len = 0;
    b->s = malloc(b->cap);
    if (b->s) b->s[0] = '\0';
}

/* Every append is a no-op once allocation has failed, so callers never have to
 * check; the failure is reported once, by returning NULL at the end. */
static bool db_ensure(DBuf* b, size_t extra) {
    if (!b->s) return false;
    if (b->len + extra + 1 <= b->cap) return true;
    size_t cap = b->cap;
    while (cap < b->len + extra + 1) cap *= 2;
    char* t = realloc(b->s, cap);
    if (!t) { free(b->s); b->s = NULL; return false; }
    b->s = t; b->cap = cap;
    return true;
}

static void db_cat(DBuf* b, const char* t) {
    size_t n = strlen(t);
    if (!db_ensure(b, n)) return;
    memcpy(b->s + b->len, t, n + 1);
    b->len += n;
}

static void db_catf(DBuf* b, const char* fmt, ...) {
    if (!b->s) return;
    va_list ap, ap2;
    va_start(ap, fmt);
    va_copy(ap2, ap);
    int need = vsnprintf(NULL, 0, fmt, ap);
    va_end(ap);
    if (need > 0 && db_ensure(b, (size_t)need)) {
        vsnprintf(b->s + b->len, (size_t)need + 1, fmt, ap2);
        b->len += (size_t)need;
    }
    va_end(ap2);
}

/* Pad with spaces to `col` counted from the last newline, so the rendered
 * column lines up however wide the operands came out. */
static void db_pad(DBuf* b, size_t col) {
    if (!b->s) return;
    size_t bol = b->len;
    while (bol > 0 && b->s[bol - 1] != '\n') bol--;
    size_t w = b->len - bol;
    while (w++ < col) db_cat(b, " ");
}

/* ------------------------------------------------------------------ *
 *  Small formatters                                                   *
 * ------------------------------------------------------------------ */

/* Shortest form that reads back exactly: 15 significant digits covers almost
 * every constant a user writes, and 17 is the guaranteed double round trip.
 *
 * A trailing "." is forced on whole numbers.  This is a TYPED machine: a Real
 * 1 and an Integer 1 live in different registers and feed different opcodes,
 * and printing both as "1" would hide exactly the coercion bugs someone reads a
 * disassembly to find. */
static void fmt_real(char* out, size_t n, double v) {
    snprintf(out, n, "%.15g", v);
    if (strtod(out, NULL) != v) snprintf(out, n, "%.17g", v);
    if (!strpbrk(out, ".eEnif")) {          /* no point, exponent, nan or inf */
        size_t len = strlen(out);
        if (len + 2 <= n) { out[len] = '.'; out[len + 1] = '\0'; }
    }
}

static void fmt_complex(char* out, size_t n, double _Complex z) {
    char re[32], im[32];       /* %.17g of a double is at most 24 characters */
    double i = cimag(z);
    fmt_real(re, sizeof re, creal(z));
    fmt_real(im, sizeof im, i < 0 ? -i : i);
    snprintf(out, n, "%s %c %s I", re, i < 0 ? '-' : '+', im);
}

static const char* ct_scalar_name(CompileType t) {
    switch (t) {
        case CT_BOOL:    return "Boolean";
        case CT_INT:     return "Integer";
        case CT_REAL:    return "Real";
        case CT_COMPLEX: return "Complex";
        default: break;
    }
    return "Unknown";
}

/* Unlike ct_name in compiled_function.c, which flattens every array to
 * "Array", this keeps the element type and rank — the two things you need to
 * read an array program's registers. */
static void ct_label(CompileType t, char* out, size_t n) {
    if (CT_IS_ARRAY(t)) snprintf(out, n, "%s[%d]", ct_scalar_name(CT_ELEM(t)), CT_RANK(t));
    else                snprintf(out, n, "%s", ct_scalar_name(t));
}

/* Named by BANK, numbered by frame slot: R scalar, V array handle, T tile.  An
 * array ARGUMENT is the one register that holds a handle while living below
 * arr_base (arguments always occupy [0, nargs)), so its declared type decides,
 * not its address — otherwise the same handle would read as `R0` here and `V4`
 * two lines later. */
static void reg_name(const CompiledProgram* p, uint32_t r, char* out, size_t n) {
    if ((int)r >= p->tile_base)     snprintf(out, n, "T%u", r);
    else if ((int)r >= p->arr_base) snprintf(out, n, "V%u", r);
    else if (r < p->nargs && p->arg_types && CT_IS_ARRAY(p->arg_types[r]))
                                    snprintf(out, n, "V%u", r);
    else                            snprintf(out, n, "R%u", r);
}

/* ------------------------------------------------------------------ *
 *  Opcode-name analysis                                               *
 * ------------------------------------------------------------------ */
/* The `_I` / `_R` / `_C` suffix of a monomorphic opcode states the type it
 * operates on: the FIRST type letter is the operand type and the LAST the
 * result, so KERN_RC is real in, complex out.  `K` marks the folded constant of
 * a K_BINK form and is not a type; digits appear in KERN_R2R and I2C.
 *
 * The tail must be SHORT and made only of those characters, or names like
 * V_KERN2 ("R" inside "KERN") and A_SHAPECHK ("C" inside "CHK") would be read
 * as typed when they are not.
 */
static CompileType suffix_type(const char* name, bool want_result) {
    const char* u = strrchr(name, '_');
    if (!u) return CT_ERR;
    const char* tail = u + 1;
    size_t len = strlen(tail);
    if (len == 0 || len > 3) return CT_ERR;
    CompileType last = CT_ERR;
    for (size_t i = 0; i < len; i++) {
        char q = tail[i];
        CompileType t = (q == 'I') ? CT_INT : (q == 'R') ? CT_REAL
                      : (q == 'C') ? CT_COMPLEX : CT_ERR;
        if (t == CT_ERR) {
            if (q == 'K' || (q >= '0' && q <= '9')) continue;
            return CT_ERR;                       /* not a type suffix at all */
        }
        if (!want_result) return t;
        last = t;
    }
    return last;
}

/* Operator base: the name before the type suffix, with the tile-form `V`
 * prefix stripped, so ADD_I / ADD_R / ADD_RK / VADD_R share one table entry.
 * `V_EW` and friends keep their `V` — the underscore says they are array ops,
 * not tile forms. */
static void op_base(const char* name, char* out, size_t n) {
    if (name[0] == 'V' && name[1] && name[1] != '_') name++;
    size_t i = 0;
    while (name[i] && name[i] != '_' && i + 1 < n) { out[i] = name[i]; i++; }
    out[i] = '\0';
}

typedef struct { const char* base; const char* glyph; } BaseText;

static const BaseText binop_glyph[] = {
    { "ADD", " + " },{ "SUB", " - " }, { "MUL", " * " }, { "DIV", " / " },
    { "POW", "^" },  { "MOD", " mod " }, { "QUOT", " quot " },
    { "LT", " < " }, { "LE", " <= " }, { "GT", " > " },  { "GE", " >= " },
    { "EQ", " == " },{ "NE", " != " }, { "AND", " && " },{ "OR", " || " },
    { "XOR", " xor " },
};

static const BaseText unop_fn[] = {
    { "SQRT", "Sqrt" },  { "EXP", "Exp" },     { "LOG", "Log" },
    { "SIN", "Sin" },    { "COS", "Cos" },     { "TAN", "Tan" },
    { "SINH", "Sinh" },  { "COSH", "Cosh" },   { "TANH", "Tanh" },
    { "ASIN", "ArcSin" },{ "ACOS", "ArcCos" }, { "ATAN", "ArcTan" },
    { "ABS", "Abs" },    { "SIGN", "Sign" },   { "FLOOR", "Floor" },
    { "CEIL", "Ceiling" },{ "ROUND", "Round" },{ "TRUNC", "IntegerPart" },
    { "RE", "Re" },      { "IM", "Im" },       { "ARG", "Arg" },
    { "CONJ", "Conjugate" }, { "ERF", "Erf" }, { "ERFC", "Erfc" },
};

static const BaseText binop_fn[] = {
    { "ATAN2", "ArcTan" }, { "MAX", "Max" }, { "MIN", "Min" },
};

static const char* lookup_base(const BaseText* tab, size_t n, const char* base) {
    for (size_t i = 0; i < n; i++) if (!strcmp(tab[i].base, base)) return tab[i].glyph;
    return NULL;
}
#define NELEM(a) (sizeof (a) / sizeof (a)[0])

/* K_BINK naming: `_KR` / `_KI` put the folded constant on the LEFT of a
 * non-commutative operator (K - R), everything else on the right. */
static bool bink_const_first(const char* name) {
    size_t n = strlen(name);
    return n >= 3 && name[n - 2] == 'K' && (name[n - 1] == 'R' || name[n - 1] == 'I');
}

static bool bink_const_is_real(const char* name) {
    size_t n = strlen(name);
    if (n < 3) return false;
    char a = name[n - 2], b = name[n - 1];
    return a == 'R' || b == 'R';       /* _RK or _KR; _IK / _KI are integers */
}

/* ------------------------------------------------------------------ *
 *  Machine-kernel reverse lookup                                      *
 * ------------------------------------------------------------------ */
/* A kernel immediate is either the descriptor registered on a SymbolDef
 * (V_KERN / V_KERN2) or one of the raw function pointers inside it (the scalar
 * KERN_* forms take k->real or k->cplx directly).  Either way the symbol table
 * is the only thing that knows the name, so scan it.
 *
 * O(table) per kernel instruction, which is fine for a debugging printer and
 * buys output with no addresses in it — the property that makes the result
 * stable enough to assert on in a test. */
typedef struct { const void* target; const char* name; } KernScan;

static void kern_visit(const char* name, SymbolDef* d, void* user) {
    KernScan* s = (KernScan*)user;
    if (s->name || !d) return;
    const NDUnaryKernel* u = (const NDUnaryKernel*)d->ndarray_unary_kernel;
    if (u && ((const void*)u == s->target
              || (const void*)u->cplx == s->target
              || (const void*)u->real == s->target)) { s->name = name; return; }
    const NDBinaryKernel* b = (const NDBinaryKernel*)d->ndarray_binary_kernel;
    if (b && ((const void*)b == s->target
              || (const void*)b->cplx == s->target)) { s->name = name; return; }
    const NDNaryKernel* k = (const NDNaryKernel*)d->ndarray_nary_kernel;
    if (k && ((const void*)k == s->target
              || (const void*)k->cplx == s->target)) { s->name = name; return; }
}

static const char* kernel_name(const void* target) {
    if (!target) return NULL;
    KernScan s; s.target = target; s.name = NULL;
    symtab_for_each(kern_visit, &s);
    return s.name;
}

/* ------------------------------------------------------------------ *
 *  Which Slot member an immediate uses                                *
 * ------------------------------------------------------------------ */
typedef enum {
    IMM_NONE, IMM_INT, IMM_REAL, IMM_KERNEL, IMM_CALL, IMM_PLOOP,
    IMM_PARTSPEC, IMM_CONST
} ImmKind;

static ImmKind imm_kind(uint16_t op) {
    switch (op) {
        case OP_CONST:                    return IMM_CONST;
        case OP_CALL:                     return IMM_CALL;
        case OP_APAR:                     return IMM_PLOOP;
        case OP_A_PART: case OP_A_PARTSET:return IMM_PARTSPEC;
        case OP_A_AXIS: case OP_A_NEW: case OP_V_EW: return IMM_INT;
        case OP_V_KERN: case OP_V_KERN2:  return IMM_KERNEL;
        default: break;
    }
    switch (compile_op_kind[op]) {
        case K_KERN1: case K_KERN2: case K_NARY: return IMM_KERNEL;
        case K_POWI:  case K_INC:   case K_LOOP: return IMM_INT;
        case K_BINK: return bink_const_is_real(op_name[op]) ? IMM_REAL : IMM_INT;
        default: break;
    }
    return IMM_NONE;
}

/* OP_CONST copies the whole Slot and carries no type tag; the register's static
 * type does not survive into the finished program.  Recover it from the first
 * instruction that READS the destination — the position it is read in says what
 * it must be.  A CONST with no reader is dead code the optimiser removes, so
 * the CT_ERR fallback (raw bits) is effectively unreachable.
 *
 * MOVE is followed rather than treated as a reader: initialising a Module local
 * puts a MOVE between the constant and every instruction that says anything
 * about its type, and `MOVE` has no type suffix of its own.
 *
 * `*dead` reports that NOTHING reads the destination before it is written
 * again — a store the optimiser left behind.  Worth distinguishing from
 * "read, but by an opcode that does not name a type": one is a missed
 * optimisation to look at, the other is a limit of this analysis. */
static CompileType const_type(const CompiledProgram* p, size_t at, bool* dead) {
    uint32_t d = p->code[at].dst;
    *dead = true;
    for (size_t j = at + 1; j < p->n; j++) {
        const Instr* c = &p->code[j];
        unsigned k = compile_op_kind[c->op];

        if (k == K_CALL || k == K_NARY || k == K_ARANGE) {
            if (d >= c->a && d < c->a + c->flags) {
                *dead = false;
                if (c->op == OP_CALL) {
                    const CompiledProgram* q = (const CompiledProgram*)c->imm.p;
                    size_t idx = d - c->a;
                    if (q && idx < q->nargs) return q->arg_types[idx];
                    return CT_ERR;
                }
                return (k == K_ARANGE) ? CT_INT : CT_ERR;   /* dims / subscripts */
            }
        }

        /* Positions that are an index or a count whatever the opcode's suffix. */
        bool idx = ((k == K_LOOP || k == K_APAR) && c->a == d)
                || (k == K_AIDX && c->a == d)
                || (k == K_ALOAD && c->b == d)
                || ((c->op == OP_A_STORE_R || c->op == OP_A_STORE_C) && c->a == d)
                || ((c->op == OP_VSETLEN || c->op == OP_VLOAD_R || c->op == OP_VLOAD_C
                     || c->op == OP_VSTORE_R || c->op == OP_VSTORE_C)
                    && (c->a == d || c->b == d));
        if (idx) { *dead = false; return CT_INT; }
        bool bl = (k == K_JZ && c->a == d)
               || ((c->op == OP_AND || c->op == OP_OR || c->op == OP_XOR
                    || c->op == OP_NOT) && (c->a == d || c->b == d));
        if (bl) { *dead = false; return CT_BOOL; }

        /* Read-modify-write positions. */
        if (c->dst == d) {
            if (k == K_RET) { *dead = false; return p->result_type; }
            if (k == K_INC || k == K_LOOP || k == K_AIDX || k == K_APAR) {
                *dead = false; return CT_INT;
            }
            if (k == K_VACC) {
                *dead = false;
                CompileType t = suffix_type(op_name[c->op], false);
                if (t != CT_ERR) return t;
            }
        }

        /* Follow the copy: whatever types the destination types this constant. */
        if (k == K_MOVE && c->a == d) { *dead = false; d = c->dst; continue; }

        /* Ordinary operand positions: the opcode's own type suffix. */
        bool reads_a = (k != K_CONST && k != K_JMP && k != K_NOP && k != K_INC
                        && k != K_RET);
        bool reads_b = (k == K_BIN || k == K_KERN2 || k == K_ALOAD || k == K_AIDX
                        || k == K_ARR || k == K_ASTORE);
        if ((reads_a && c->a == d) || (reads_b && c->b == d)) {
            *dead = false;
            CompileType t = suffix_type(op_name[c->op], false);
            if (t != CT_ERR) return t;
        }

        /* Redefined before anything read it — stop, or a later unrelated use of
         * the same register would be mistaken for this constant's. */
        if (c->dst == d && k != K_JMP && k != K_JZ && k != K_RET
            && k != K_ASTORE && k != K_NOP) break;
    }
    return CT_ERR;
}

/* ------------------------------------------------------------------ *
 *  Callee-program worklist                                           *
 * ------------------------------------------------------------------ */
/* OP_CALL's immediate is a BORROWED callee program.  Numbering them and
 * listing each once is also what makes mutual recursion (A compiled against B,
 * B later compiled against A) terminate. */
typedef struct { const CompiledProgram** v; int n, cap; } ProgList;

static int pl_index(ProgList* L, const CompiledProgram* q) {
    for (int i = 0; i < L->n; i++) if (L->v[i] == q) return i;
    if (L->n == L->cap) {
        int cap = L->cap ? L->cap * 2 : 4;
        const CompiledProgram** v = realloc(L->v, (size_t)cap * sizeof *v);
        if (!v) return -1;
        L->v = v; L->cap = cap;
    }
    L->v[L->n] = q;
    return L->n++;
}

/* ------------------------------------------------------------------ *
 *  Immediate rendering                                                *
 * ------------------------------------------------------------------ */

static void render_partspec(DBuf* b, const CompiledProgram* p, const PartSpec* ps) {
    if (!ps) { db_cat(b, "<spec>"); return; }
    db_cat(b, "[[");
    for (int i = 0; i < ps->n; i++) {
        if (i) db_cat(b, ", ");
        if (ps->lit && ps->lit[i]) {
            char* s = expr_to_string(ps->lit[i]);
            db_cat(b, s ? s : "?");
            free(s);
        } else if (ps->reg && ps->reg[i] >= 0) {
            char r[24]; reg_name(p, (uint32_t)ps->reg[i], r, sizeof r);
            db_cat(b, r);
        } else {
            db_cat(b, "?");
        }
    }
    db_cat(b, "]]");
}

/* Writes the immediate of instruction `i`.  Returns false when the opcode has
 * no immediate, so the caller can skip the separator. */
static bool render_imm(DBuf* b, const CompiledProgram* p, size_t i, ProgList* L) {
    const Instr* c = &p->code[i];
    char tmp[80];
    switch (imm_kind(c->op)) {
        case IMM_NONE: return false;
        case IMM_INT:  db_catf(b, "%lld", c->imm.i); return true;
        case IMM_REAL: fmt_real(tmp, sizeof tmp, c->imm.r); db_cat(b, tmp); return true;
        case IMM_KERNEL: {
            const char* nm = kernel_name(c->imm.p);
            db_catf(b, "<%s>", nm ? nm : "kernel");
            return true;
        }
        case IMM_CALL: {
            int k = pl_index(L, (const CompiledProgram*)c->imm.p);
            if (k < 0) db_cat(b, "<call>"); else db_catf(b, "<call #%d>", k);
            return true;
        }
        case IMM_PLOOP: {
            const ParLoop* pl = (const ParLoop*)c->imm.p;
            if (pl && p->ploops) db_catf(b, "<ploop #%d>", (int)(pl - p->ploops));
            else                 db_cat(b, "<ploop>");
            return true;
        }
        case IMM_PARTSPEC:
            render_partspec(b, p, (const PartSpec*)c->imm.p);
            return true;
        case IMM_CONST: {
            bool dead = false;
            switch (const_type(p, i, &dead)) {
                case CT_BOOL:    db_cat(b, c->imm.i ? "True" : "False"); break;
                case CT_INT:     db_catf(b, "%lld", c->imm.i); break;
                case CT_REAL:    fmt_real(tmp, sizeof tmp, c->imm.r); db_cat(b, tmp); break;
                case CT_COMPLEX: fmt_complex(tmp, sizeof tmp, c->imm.z); db_cat(b, tmp); break;
                /* Nothing reads it before it is written again: the value has no
                 * type because it has no meaning.  Saying so is more use than
                 * the bit pattern — it names a store the optimiser left in. */
                default:         if (dead) db_cat(b, "<dead store>");
                                 else db_catf(b, "0x%016llx", (unsigned long long)c->imm.i);
                                 break;
            }
            return true;
        }
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Operand column                                                     *
 * ------------------------------------------------------------------ */
/* `b` is a branch TARGET rather than a register on exactly these kinds — the
 * same predicate patch_reg uses, asked of the same table, so the two cannot
 * disagree about what a field means. */
static bool op_b_is_target(uint16_t op) {
    unsigned k = compile_op_kind[op];
    return k == K_JMP || k == K_JZ || k == K_LOOP || k == K_APAR;
}

static void render_operands(DBuf* b, const CompiledProgram* p, size_t i, ProgList* L) {
    const Instr* c = &p->code[i];
    unsigned k = compile_op_kind[c->op];
    char rd[24], ra[24], rb[24];
    reg_name(p, c->dst, rd, sizeof rd);
    reg_name(p, c->a, ra, sizeof ra);
    reg_name(p, c->b, rb, sizeof rb);

    switch (k) {
        case K_NOP:
            break;
        case K_CONST:
            db_catf(b, "%s, ", rd);
            render_imm(b, p, i, L);
            break;
        case K_JMP:
            db_catf(b, "-> %u", c->b);
            break;
        case K_JZ:
            db_catf(b, "%s, -> %u", ra, c->b);
            break;
        case K_LOOP: case K_APAR:
            db_catf(b, "%s, %s, -> %u", rd, ra, c->b);
            if (k == K_APAR) { db_cat(b, "  "); render_imm(b, p, i, L); }
            break;
        case K_RET:
            db_cat(b, rd);
            break;
        case K_INC:
            db_catf(b, "%s, ", rd);
            render_imm(b, p, i, L);
            break;
        case K_CALL: case K_NARY: case K_ARANGE: {
            db_catf(b, "%s, [", rd);
            for (unsigned j = 0; j < c->flags; j++) {
                char r[24]; reg_name(p, c->a + j, r, sizeof r);
                db_catf(b, "%s%s", j ? ", " : "", r);
            }
            db_cat(b, "]");
            /* A_PART reads the source array through `b`, outside the range. */
            if (c->op == OP_A_PART) db_catf(b, ", %s", rb);
            if (c->op == OP_A_PARTSET) db_catf(b, ", %s", rb);
            db_cat(b, "  ");
            render_imm(b, p, i, L);
            break;
        }
        default: {
            /* Fixed-field forms: dst, a, and b when it is a register.  The two
             * pure checks write nothing, so their `dst` field is never set and
             * printing it would name an unrelated register. */
            bool has_dst = (c->op != OP_VSETLEN && c->op != OP_A_SHAPECHK);
            if (has_dst) db_catf(b, "%s, ", rd);
            db_cat(b, ra);
            bool has_b = (k == K_BIN || k == K_KERN2 || k == K_ALOAD || k == K_AIDX
                          || k == K_ARR || k == K_ASTORE);
            if (has_b && !op_b_is_target(c->op)) db_catf(b, ", %s", rb);
            if (imm_kind(c->op) != IMM_NONE) { db_cat(b, ", "); render_imm(b, p, i, L); }
            /* K_ARR is the only kind whose flags carry ownership and operand
             * kinds; for K_ASTORE and the range kinds they mean something else
             * and are shown elsewhere. */
            if (k == K_ARR && c->flags) {
                static const char* akname[] = { "arr", "real", "cplx" };
                unsigned f = c->flags, ka = AF_A(f), kb = AF_B(f);
                db_cat(b, "  [");
                if (f & AF_FREE_A) db_cat(b, "free a; ");
                if (f & AF_FREE_B) db_cat(b, "free b; ");
                db_catf(b, "a:%s b:%s -> %s",
                        ka < 3 ? akname[ka] : "?", kb < 3 ? akname[kb] : "?",
                        ct_scalar_name((CompileType)AF_R(f)));
                db_cat(b, "]");
            }
            break;
        }
    }
}

/* ------------------------------------------------------------------ *
 *  Rendered (pseudo-source) column                                    *
 * ------------------------------------------------------------------ */
static void render_meaning(DBuf* b, const CompiledProgram* p, size_t i, ProgList* L) {
    const Instr* c = &p->code[i];
    unsigned k = compile_op_kind[c->op];
    const char* nm = op_name[c->op];
    char base[32]; op_base(nm, base, sizeof base);
    char rd[24], ra[24], rb[24];
    reg_name(p, c->dst, rd, sizeof rd);
    reg_name(p, c->a, ra, sizeof ra);
    reg_name(p, c->b, rb, sizeof rb);

    switch (k) {
        case K_NOP:   db_cat(b, "(removed)"); return;
        case K_CONST: db_catf(b, "%s = ", rd); render_imm(b, p, i, L); return;
        case K_MOVE:  db_catf(b, "%s = %s", rd, ra); return;
        case K_JMP:   db_catf(b, "goto %u", c->b); return;
        case K_JZ:    db_catf(b, "if !%s goto %u", ra, c->b); return;
        case K_LOOP:  db_catf(b, "if ++%s < %s goto %u", rd, ra, c->b); return;
        case K_APAR:  db_catf(b, "parallel %s over [0, %s), else fall through; join %u",
                              rd, ra, c->b); return;
        case K_RET:   db_catf(b, "return %s", rd); return;
        case K_INC:   db_catf(b, "%s = %s + %lld", rd, rd, c->imm.i); return;
        case K_POWI:  db_catf(b, "%s = %s^%lld", rd, ra, c->imm.i); return;
        case K_AIDX:  db_catf(b, "%s = %s*dim(%s, %lld) + resolve(%s)",
                              rd, rd, rb, c->imm.i, ra); return;
        case K_ALOAD: db_catf(b, "%s = %s[%s]", rd, ra, rb); return;
        case K_VACC:  db_catf(b, "%s += total(%s)", rd, ra); return;

        case K_KERN1: case K_KERN2: case K_NARY: {
            const char* fn = kernel_name(c->imm.p);
            db_catf(b, "%s = %s[", rd, fn ? fn : "kernel");
            if (k == K_NARY) {
                for (unsigned j = 0; j < c->flags; j++) {
                    char r[24]; reg_name(p, c->a + j, r, sizeof r);
                    db_catf(b, "%s%s", j ? ", " : "", r);
                }
            } else if (k == K_KERN2) {
                db_catf(b, "%s, %s", ra, rb);
            } else {
                db_cat(b, ra);
            }
            db_cat(b, "]");
            return;
        }

        case K_CALL: {
            int idx = pl_index(L, (const CompiledProgram*)c->imm.p);
            db_catf(b, "%s = program#%d(", rd, idx);
            for (unsigned j = 0; j < c->flags; j++) {
                char r[24]; reg_name(p, c->a + j, r, sizeof r);
                db_catf(b, "%s%s", j ? ", " : "", r);
            }
            db_cat(b, ")");
            return;
        }

        case K_BINK: {
            const char* g = lookup_base(binop_glyph, NELEM(binop_glyph), base);
            char kv[80];
            if (imm_kind(c->op) == IMM_REAL) fmt_real(kv, sizeof kv, c->imm.r);
            else snprintf(kv, sizeof kv, "%lld", c->imm.i);
            if (g) {
                if (bink_const_first(nm)) db_catf(b, "%s = %s%s%s", rd, kv, g, ra);
                else                      db_catf(b, "%s = %s%s%s", rd, ra, g, kv);
            } else {
                db_catf(b, "%s = %s(%s, %s)", rd, nm, ra, kv);
            }
            return;
        }

        case K_BIN: {
            const char* g = lookup_base(binop_glyph, NELEM(binop_glyph), base);
            if (g) { db_catf(b, "%s = %s%s%s", rd, ra, g, rb); return; }
            const char* fn = lookup_base(binop_fn, NELEM(binop_fn), base);
            if (fn) { db_catf(b, "%s = %s[%s, %s]", rd, fn, ra, rb); return; }
            if (!strcmp(base, "LOAD")) { db_catf(b, "%s = %s[%s ...]", rd, ra, rb); return; }
            db_catf(b, "%s = %s(%s, %s)", rd, nm, ra, rb);
            return;
        }

        case K_UN: {
            if (!strcmp(base, "NEG")) { db_catf(b, "%s = -%s", rd, ra); return; }
            if (!strcmp(base, "INV")) { db_catf(b, "%s = 1/%s", rd, ra); return; }
            if (!strcmp(base, "NOT")) { db_catf(b, "%s = !%s", rd, ra); return; }
            if (c->op == OP_A_SIZE) { db_catf(b, "%s = Length[%s] (flat)", rd, ra); return; }
            if (!strcmp(base, "I2R") || !strcmp(base, "I2C") || !strcmp(base, "R2C")) {
                /* widening coercion: the last letter of the base names the target */
                char to = base[2];
                db_catf(b, "%s = (%s) %s", rd,
                        ct_scalar_name(to == 'C' ? CT_COMPLEX : CT_REAL), ra);
                return;
            }
            if (!strcmp(base, "SPLAT")) { db_catf(b, "%s = splat(%s)", rd, ra); return; }
            const char* fn = lookup_base(unop_fn, NELEM(unop_fn), base);
            if (fn) { db_catf(b, "%s = %s[%s]", rd, fn, ra); return; }
            db_catf(b, "%s = %s(%s)", rd, nm, ra);
            return;
        }

        case K_ARR: case K_ASTORE: {
            switch (c->op) {
                case OP_ARR_FREE:  db_catf(b, "free %s", rd); return;
                case OP_A_COPY:    db_catf(b, "%s = copy(%s)", rd, ra); return;
                case OP_A_XFER:    db_catf(b, "%s = move(%s)", rd, ra); return;
                case OP_A_NEWLIKE: db_catf(b, "%s = buffer like %s (%s)", rd, ra,
                                           ct_scalar_name((CompileType)AF_R(c->flags))); return;
                case OP_A_SHAPECHK:db_catf(b, "assert dims(%s) == dims(%s)", ra, rb); return;
                case OP_V_LEN:     db_catf(b, "%s = Length[%s]", rd, ra); return;
                case OP_V_TOTAL:   db_catf(b, "%s = Total[%s]", rd, ra); return;
                case OP_V_EW:      db_catf(b, "%s = %s %s %s", rd, ra,
                                           c->imm.i ? "+" : "*", rb); return;
                case OP_V_POW:     db_catf(b, "%s = %s^%s", rd, ra, rb); return;
                case OP_V_KERN: case OP_V_KERN2: {
                    const char* fn = kernel_name(c->imm.p);
                    if (c->op == OP_V_KERN) db_catf(b, "%s = %s[%s]", rd, fn ? fn : "kernel", ra);
                    else db_catf(b, "%s = %s[%s, %s]", rd, fn ? fn : "kernel", ra, rb);
                    return;
                }
                case OP_A_STORE_R: case OP_A_STORE_C:
                    db_catf(b, "%s[%s] = %s", rd, ra, rb); return;
                case OP_VSETLEN:
                    db_catf(b, "vlen = min(%s - %s, %d)", rb, ra, VBLOCK); return;
                case OP_VSTORE_R: case OP_VSTORE_C:
                    db_catf(b, "%s[%s ...] = %s", rd, ra, rb); return;
                default: break;
            }
            db_catf(b, "%s = %s(%s, %s)", rd, nm, ra, rb);
            return;
        }

        case K_ARANGE: {
            if (c->op == OP_A_NEW) {
                db_catf(b, "%s = zeros(", rd);
                for (unsigned j = 0; j < c->flags; j++) {
                    char r[24]; reg_name(p, c->a + j, r, sizeof r);
                    db_catf(b, "%s%s", j ? " x " : "", r);
                }
                db_catf(b, ") of %s", ct_scalar_name((CompileType)c->imm.i));
                return;
            }
            if (c->op == OP_A_PART) {
                db_catf(b, "%s = %s", rd, rb);
                render_partspec(b, p, (const PartSpec*)c->imm.p);
                return;
            }
            if (c->op == OP_A_PARTSET) {
                db_cat(b, rd);
                render_partspec(b, p, (const PartSpec*)c->imm.p);
                db_catf(b, " = %s", rb);
                return;
            }
            db_cat(b, nm);
            return;
        }
    }
    db_cat(b, nm);
}

/* ------------------------------------------------------------------ *
 *  Program listing                                                    *
 * ------------------------------------------------------------------ */
#define MEANING_COL 48

static void emit_listing(DBuf* b, const CompiledProgram* p, ProgList* L) {
    /* Mark branch targets so loop structure is visible without chasing numbers. */
    unsigned char* tgt = calloc(p->n ? p->n : 1, 1);
    if (tgt)
        for (size_t i = 0; i < p->n; i++)
            if (op_b_is_target(p->code[i].op) && p->code[i].b < p->n)
                tgt[p->code[i].b] = 1;

    for (size_t i = 0; i < p->n; i++) {
        db_catf(b, "%s%4zu  %-10s ", (tgt && tgt[i]) ? ">" : " ", i, op_name[p->code[i].op]);
        render_operands(b, p, i, L);
        db_pad(b, MEANING_COL);
        db_cat(b, "  ");
        render_meaning(b, p, i, L);
        db_cat(b, "\n");
    }
    free(tgt);
}

static void emit_header(DBuf* b, const CompiledProgram* p, const char* const* names) {
    char ty[48], rr[24];
    if (p->nargs == 0) {
        db_cat(b, "Arguments   none\n");
    } else {
        db_catf(b, "Arguments   %zu\n", p->nargs);
        for (size_t i = 0; i < p->nargs; i++) {
            char r[24];
            reg_name(p, (uint32_t)i, r, sizeof r);
            ct_label(p->arg_types[i], ty, sizeof ty);
            db_catf(b, "              %-4s : %-12s", r, ty);
            if (names && names[i]) db_catf(b, " %s", names[i]);
            /* argdep is what the colored-FD Jacobian reads; an argument the body
             * never touches is worth seeing here too. */
            if (p->argdep && !p->argdep[i]) db_cat(b, "   (unused)");
            db_cat(b, "\n");
        }
    }
    reg_name(p, (uint32_t)p->result_reg, rr, sizeof rr);
    ct_label(p->result_type, ty, sizeof ty);
    db_catf(b, "Result      %s : %s\n", rr, ty);
    db_catf(b, "Registers   %d scalar, %d array, %d tile   (frame %zu slots)\n",
            p->arr_base, p->tile_base - p->arr_base, p->nreg - p->tile_base,
            p->frame_slots);
    db_catf(b, "Program     %zu instruction%s, %d CSE", p->n, p->n == 1 ? "" : "s", p->ncse);
    if (p->nploops) db_catf(b, ", %d parallel loop%s", p->nploops,
                            p->nploops == 1 ? "" : "s");
    if (p->all_real) db_cat(b, ", all-Real fast path");
    db_cat(b, "\n\n");
}

char* compiled_disassemble(const CompiledProgram* p, const char* const* arg_names) {
    if (!p) return NULL;
    DBuf b; db_init(&b);
    ProgList L; L.v = NULL; L.n = 0; L.cap = 0;
    if (pl_index(&L, p) < 0) { free(b.s); return NULL; }

    /* L grows as OP_CALLs are discovered, so this is a worklist: every callee
     * reached from the entry program is listed exactly once. */
    for (int i = 0; i < L.n; i++) {
        if (i) db_catf(&b, "\n--- called program #%d ---\n\n", i);
        emit_header(&b, L.v[i], i ? NULL : arg_names);
        emit_listing(&b, L.v[i], &L);
    }
    free(L.v);
    return b.s;
}

/* ------------------------------------------------------------------ *
 *  CompiledFunction level                                             *
 * ------------------------------------------------------------------ */
char* compiled_function_disassemble(const CompiledFunction* cf) {
    if (!cf) return NULL;
    DBuf b; db_init(&b);

    size_t n = compiled_function_num_args(cf);
    const char* const* names = compiled_function_arg_names(cf);
    const CompileType* types = compiled_function_arg_types(cf);
    const Expr* body = compiled_function_body(cf);

    db_cat(&b, "Signature   CompiledFunction[{");
    for (size_t i = 0; i < n; i++) {
        char ty[48];
        ct_label(types ? types[i] : CT_REAL, ty, sizeof ty);
        db_catf(&b, "%s%s : %s", i ? ", " : "",
                (names && names[i]) ? names[i] : "?", ty);
    }
    db_cat(&b, "}, ");
    char* bs = expr_to_string((Expr*)body);
    db_cat(&b, bs ? bs : "?");
    free(bs);
    db_cat(&b, "]\n");

    const CompiledProgram* p = compiled_function_program(cf);
    if (p) {
        char* prog = compiled_disassemble(p, names);
        if (prog) { db_cat(&b, prog); free(prog); }
        else if (b.s) { free(b.s); b.s = NULL; }
    } else {
        /* Not compiled: the useful answer is not "no bytecode" but WHY, and the
         * bail globals are only valid until the next compile call — so re-run
         * the compile to repopulate them, exactly as CompileDiagnostics does. */
        db_cat(&b, "Program     not compiled — every call runs the interpreter\n");
        CompiledProgram* again = compile_expr(body, names, types, n);
        if (again) {
            /* Compiles now but did not at construction: report rather than
             * silently print a program the object will never actually run. */
            db_cat(&b, "            (the body compiles now; this object predates that)\n");
            compiled_free(again);
        } else {
            const char* why  = compiled_bail_reason();
            const char* what = compiled_bail_expr();
            db_catf(&b, "Reason      %s\n", why ? why : "unknown");
            if (what) db_catf(&b, "Bailed on   %s\n", what);
        }
    }
    return b.s;
}
