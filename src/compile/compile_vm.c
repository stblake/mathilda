/* Mathilda — Compile[]: the bytecode virtual machine.
 *
 * Executes a finished CompiledProgram over a flat Slot register file.  The
 * instruction set is monomorphic — the opcode carries the operand type, so the
 * hot loop does no tag dispatch and the Real path never pays complex cost.  This
 * file holds the scalar integer/boxing helpers, the association and array
 * runtime, the core vm_run / vm_call / parallel-strip machinery, the
 * arbitrary-precision managed arena, the per-call runtime frame, and the public
 * eval entry points (compiled_eval / _real / _batch, compiled_free) with the
 * program accessors.
 *
 * It shares nothing with the emitter but the finished-program representation in
 * compile_internal.h — it never sees the Ctx builder or compile_emit.h — plus a
 * handful of leaf helpers the emit modules also use (the assoc value coercions
 * and the managed-constant materialisers), declared in compile_internal.h. */
#include "compile.h"
#include "compile_internal.h"    /* Slot / Instr / CompiledProgram / opcodes / VM_TLS */
#include "../arithmetic.h"       /* expr_new_* / expr_free / expr_bigint_normalize */
#include "../symtab.h"
#include "../ndarray.h"          /* ndarray_part(_set) / map / elementwise / NDType */
#include "../assoc.h"            /* assoc_lookup_value / assoc_values_list / assoc_set_key ... */
#include "../ndreduce.h"         /* ndred_total_all */
#include "../ndarray_internal.h" /* nd_parallel_for — threading the fused map loop */
#include "../sym_names.h"        /* SYM_Association / SYM_All / SYM_Span / SYM_List */
#ifdef USE_MPFR
#include "../numeric_complex.h"  /* ncpx / numeric_mpfr_* — the managed containers */
#endif
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.7182818284590452354
#endif

/* The integer form is `ci_powi` in compile_internal.h: same binary exponentiation,
 * overflow-checked at every step, and shared with the optimiser's folding. */
static double    ipow_r(double b, long long n) { if (n < 0) { b = 1.0 / b; n = -n; } double r = 1; while (n) { if (n & 1) r *= b; b *= b; n >>= 1; } return r; }

/* ------------------------------------------------------------------ *
 *  Exact integer kernels for the integer-closed heads                 *
 * ------------------------------------------------------------------ *
 * Each returns TRUE ON FAILURE — overflow, or an argument outside the domain
 * where the interpreter yields an integer — matching the ci_* convention, so a
 * VM body is the same one-line IOP() as an arithmetic opcode.  Failure hands the
 * call to the interpreter, which then answers exactly (a bigint) or symbolically
 * (ComplexInfinity, a Rational), whichever is right.
 *
 * The domains are not guesses; they are what the interpreter actually does:
 *   Factorial[-1]     ComplexInfinity     Gamma[0]        ComplexInfinity
 *   Pochhammer[7,-3]  1/120  (Rational)   7^-3            1/343  (Rational)
 *   Arg[-3]           Pi     (Symbol)     Binomial[3,7]   0      (all integer)
 * so Factorial needs n >= 0, Gamma n >= 1, Pochhammer n >= 0, Power e >= 0, and
 * Binomial no guard at all. */

static bool int_factorial(long long n, long long* out) {
    if (n < 0) return true;                       /* ComplexInfinity */
    long long r = 1;
    for (long long k = 2; k <= n; k++) if (ci_mul(r, k, &r)) return true;
    *out = r;
    return false;                                 /* 21! already overflows */
}

/* Gamma[n] = (n-1)! on the positive integers; Gamma[0] and Gamma of a negative
 * integer are ComplexInfinity, which is not a machine number. */
static bool int_gamma(long long n, long long* out) {
    if (n < 1) return true;
    return int_factorial(n - 1, out);
}

/* Binomial[n, k] for machine integers, by the multiplicative recurrence
 * C(n,k) = C(n,k-1) * (n-k+1) / k, which is exact at every step because the
 * partial product of k consecutive integers is divisible by k!.  Computing
 * n!/(k!(n-k)!) directly would overflow at n = 21 for results that fit easily. */
static bool int_binomial(long long n, long long k, long long* out) {
    if (n >= 0) {
        if (k < 0 || k > n) { *out = 0; return false; }       /* Binomial[3,7] = 0 */
        if (k > n - k) k = n - k;                             /* symmetry: fewer steps */
        long long r = 1;
        for (long long i = 1; i <= k; i++) {
            if (ci_mul(r, n - k + i, &r)) return true;
            r /= i;                                           /* exact by construction */
        }
        *out = r;
        return false;
    }
    /* Negative upper index: Binomial[-n, k] = (-1)^k Binomial[n+k-1, k].  Left to
     * the interpreter — it is rare, and the identity is easy to get subtly wrong
     * in a way no test here would catch. */
    return true;
}

/* Pochhammer[a, n] = a (a+1) ... (a+n-1); a falling/negative n gives a Rational. */
static bool int_pochhammer(long long a, long long n, long long* out) {
    if (n < 0) return true;
    long long r = 1;
    for (long long i = 0; i < n; i++) {
        long long t;
        if (ci_add(a, i, &t) || ci_mul(r, t, &r)) return true;
    }
    *out = r;
    return false;
}

/* Fibonacci / LucasL by iteration, including the negative index (the interpreter
 * defines both there: F[-n] = (-1)^(n+1) F[n], L[-n] = (-1)^n L[n]).  Iterative
 * rather than the closed form because these must be EXACT — a double loses
 * Fibonacci exactly where it starts to matter, and F[92] is the last one that
 * fits in an int64 anyway. */
static bool int_fib2(long long n, bool lucas, long long* out) {
    if (n == LLONG_MIN) return true;              /* -n would overflow; test FIRST */
    bool neg = n < 0;
    long long m = neg ? -n : n;
    long long a = lucas ? 2 : 0, b = 1;           /* (L0,L1) = (2,1); (F0,F1) = (0,1) */
    for (long long i = 0; i < m; i++) {
        long long t;
        if (ci_add(a, b, &t)) return true;
        a = b; b = t;
    }
    /* F[-n] flips sign for even n, L[-n] for odd n. */
    if (neg && ((m % 2 == 0) == !lucas)) { if (ci_neg(a, &a)) return true; }
    *out = a;
    return false;
}
static bool int_fib(long long n, long long* out)   { return int_fib2(n, false, out); }
static bool int_lucas(long long n, long long* out) { return int_fib2(n, true,  out); }

/* a^e on the integers.  A negative exponent is a Rational in the interpreter
 * (7^-3 is 1/343) and 0^0 is Indeterminate, so both abandon the call. */
static bool int_pow(long long a, long long e, long long* out) {
    if (e < 0) return true;
    if (e == 0 && a == 0) return true;            /* Indeterminate */
    return ci_powi(a, e, out);
}

/* GCD is non-negative and GCD[0, 0] is 0, matching the interpreter.  The only
 * failure is INT64_MIN, whose magnitude is not representable. */
static bool int_gcd(long long a, long long b, long long* out) {
    if (a == LLONG_MIN || b == LLONG_MIN) return true;
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { long long t = a % b; a = b; b = t; }
    *out = a;
    return false;
}

/* LCM[a, b] = |a b| / gcd, non-negative, and 0 when either side is 0.  Divide
 * BEFORE multiplying so the intermediate stays in range whenever the answer
 * does. */
static bool int_lcm(long long a, long long b, long long* out) {
    if (a == 0 || b == 0) { *out = 0; return false; }
    long long g;
    if (int_gcd(a, b, &g) || g == 0) return true;
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    return ci_mul(a / g, b, out);
}

/* Number of digits of n in the given base; IntegerLength[0] is 0. */
static bool int_ilen(long long n, long long base, long long* out) {
    if (base < 2) return true;
    if (n == LLONG_MIN) return true;
    if (n < 0) n = -n;
    long long k = 0;
    while (n > 0) { n /= base; k++; }
    *out = k;
    return false;
}

/* Largest e with base^e dividing n.  IntegerExponent[0, b] is Infinity in the
 * interpreter, which is not a machine integer. */
static bool int_iexp(long long n, long long base, long long* out) {
    if (base < 2 || n == 0) return true;
    if (n == LLONG_MIN) return true;
    if (n < 0) n = -n;
    long long k = 0;
    while (n % base == 0) { n /= base; k++; }
    *out = k;
    return false;
}

/* PowerMod[a, e, m], including the negative exponent (a modular INVERSE, which
 * exists only when gcd(a, m) is 1 — the interpreter reports that, so we defer).
 *
 * The modulus is capped at sqrt(INT64_MAX) so every intermediate product fits an
 * int64.  Going wider needs a 128-bit multiply, and `__int128` is a GNU
 * extension this file cannot use (CLAUDE.md: strict C99); a larger modulus
 * therefore hands the call to the interpreter, which has GMP. */
#define POWMOD_MAX_M 3037000499LL             /* floor(sqrt(2^63 - 1)) */
static bool int_powmod(long long a, long long e, long long m, long long* out) {
    if (m == 0 || m > POWMOD_MAX_M || m < -POWMOD_MAX_M) return true;
    if (m < 0) return true;                   /* interpreter's sign convention */
    if (m == 1) { *out = 0; return false; }
    a %= m; if (a < 0) a += m;
    if (e < 0) {
        /* Modular inverse by the extended Euclid, then the positive power. */
        long long r0 = m, r1 = a, s0 = 0, s1 = 1;
        while (r1) {
            long long q = r0 / r1;
            long long t = r0 - q * r1; r0 = r1; r1 = t;
            t = s0 - q * s1;           s0 = s1; s1 = t;
        }
        if (r0 != 1) return true;             /* not invertible */
        a = s0 % m; if (a < 0) a += m;
        if (e == LLONG_MIN) return true;      /* -e would overflow */
        e = -e;
    }
    long long r = 1;
    while (e > 0) {
        if (e & 1) r = (r * a) % m;
        e >>= 1;
        if (e > 0) a = (a * a) % m;
    }
    *out = r;
    return false;
}
static double _Complex ipow_c(double _Complex b, long long n) { if (n < 0) { b = 1.0 / b; n = -n; } double _Complex r = 1; while (n) { if (n & 1) r *= b; b *= b; n >>= 1; } return r; }

/* ------------------------------------------------------------------ *
 *  Array opcodes (M3a): delegation to the NDArray subsystem           *
 * ------------------------------------------------------------------ */

/* Read a scalar register as a (re, im) pair, per its compile-time operand kind. */
static void vm_scalar_pair(const Slot* s, unsigned kind, double* re, double* im) {
    if      (kind == AK_COMPLEX) { *re = creal(s->z);   *im = cimag(s->z); }
    else if (kind == AK_INT)     { *re = (double)s->i;  *im = 0.0; }
    else                         { *re = s->r;          *im = 0.0; }
}

/* Box a scalar register as a temporary numeric Expr, because the ND helpers
 * take Expr operands.  Two allocations amortised over a whole-buffer pass is
 * nothing, and it keeps one implementation of the broadcast and dtype-promotion
 * rules instead of a second copy here. */
static Expr* vm_box_scalar(const Slot* s, unsigned kind) {
    /* An integer register boxes as an Integer, so the ND layer sees an EXACT
     * scalar and keeps an int64 buffer exact instead of promoting it. Reading it
     * back through the double pair would defeat the whole point of AK_INT. */
    if (kind == AK_INT) return expr_new_integer((int64_t)s->i);
    double re, im;
    vm_scalar_pair(s, kind, &re, &im);
    if (im == 0.0) return expr_new_real(re);
    return make_complex(expr_new_real(re), expr_new_real(im));
}

/* Concrete (re, im) of a scalar Expr produced by an ND reduction. */
static bool vm_num_pair(const Expr* e, double* re, double* im) {
    if (!e) return false;
    if (e->type == EXPR_REAL)    { *re = e->data.real;             *im = 0.0; return true; }
    if (e->type == EXPR_INTEGER) { *re = (double)e->data.integer;  *im = 0.0; return true; }
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol.name, "Complex") == 0) {
        double ar, ai, br, bi;
        if (vm_num_pair(e->data.function.args[0], &ar, &ai)
            && vm_num_pair(e->data.function.args[1], &br, &bi)) {
            *re = ar; *im = br; return true;
        }
    }
    return false;
}

/* Write a reduction's scalar result into a register, honouring the element type
 * the program promised: a value that left the real axis where the program
 * promised a real one fails the call. */
static bool vm_write_scalar(const Expr* e, unsigned relem, Slot* d) {
    double re, im;
    if (!vm_num_pair(e, &re, &im)) return false;
    if (relem == (unsigned)CT_COMPLEX) { d->z = re + im * I; return true; }
    if (im != 0.0) return false;
    if (relem == (unsigned)CT_INT) d->i = (long long)re; else d->r = re;
    return true;
}

/* Execute one array opcode.  Returns false to abort the whole program — a shape
 * mismatch, an allocation failure, a kernel that declined an element, or a
 * result that left the promised element type — after which the caller falls
 * back to the interpreter.  Operands flagged AF_FREE_* are consumed here, AFTER
 * the op has read them, so the result may legitimately reuse an operand's
 * register; each freed slot is NULLed so an abort can never double-free. */
/* An independent copy of `x` in the program's CANONICAL dtype for `relem`.
 *
 * Every array a program owns is float64 or complex64, never float32, and that
 * is load-bearing rather than tidy: A_STORE writes the buffer at its declared
 * width, so one float32 array reaching a store would write doubles into half-
 * sized slots.  Arguments may be any dtype — they are only ever read — so the
 * narrowing is done here, at the point ownership begins. */
/* The buffer dtype a program uses for a given element type.  A program owns only
 * these four — float32 never appears, see nd_own_copy. */
static NDType ct_elem_ndt(unsigned relem) {
    return relem == (unsigned)CT_COMPLEX ? NDT_COMPLEX64
         : relem == (unsigned)CT_INT     ? NDT_INT64
         : relem == (unsigned)CT_BOOL    ? NDT_BOOL
                                         : NDT_FLOAT64;
}

static Expr* nd_own_copy(const Expr* x, unsigned relem) {
    if (!x || x->type != EXPR_NDARRAY) return NULL;
    NDType dt = ct_elem_ndt(relem);
    NDType sdt = x->data.ndarray.dtype;
    /* An integer-typed program will not silently take a float buffer: the two
     * are different element types to the interpreter as well (`Total` of a
     * float64 NDArray is a Real), so the call goes back rather than rounding.
     * Bool is the same kind of distinct type — never coerced to/from a number. */
    if ((dt == NDT_INT64) != (sdt == NDT_INT64)) return NULL;
    if ((dt == NDT_BOOL)  != (sdt == NDT_BOOL))  return NULL;
    size_t n = ndarray_size(x), esz = ndt_elem_size(dt);
    void* buf = malloc(esz * (n ? n : 1));
    if (!buf) return NULL;
    if (sdt == dt) memcpy(buf, x->data.ndarray.data, esz * n);
    else for (size_t k = 0; k < n; k++) {
        double re, im;
        ndt_get(x->data.ndarray.data, k, sdt, &re, &im);
        /* A complex source into a real program is the array form of the scalar
         * contract: the program promised real, so it fails rather than truncate. */
        if (im != 0.0 && dt != NDT_COMPLEX64) { free(buf); return NULL; }
        ndt_set(buf, k, dt, re, im);
    }
    Expr* nw = expr_new_ndarray_like(x, x->data.ndarray.rank, x->data.ndarray.dims, buf, dt);
    if (!nw) free(buf);
    return nw;
}

/* Materialise a general Part's subscript list for this call: literal specs are
 * borrowed straight from the PartSpec, computed ones are boxed from registers. */
static bool ps_build_indices(const PartSpec* ps, const Slot* R, Expr** idx) {
    for (int i = 0; i < ps->n; i++) {
        if (ps->lit[i]) { idx[i] = ps->lit[i]; continue; }
        idx[i] = expr_new_integer((int64_t)R[ps->reg[i]].i);
        if (!idx[i]) {
            for (int j = 0; j < i; j++) if (!ps->lit[j]) expr_free(idx[j]);
            return false;
        }
    }
    return true;
}
static void ps_free_indices(const PartSpec* ps, Expr** idx) {
    for (int i = 0; i < ps->n; i++) if (!ps->lit[i]) expr_free(idx[i]);
}

/* The array opcodes whose operands are a RUN of registers, so they need the
 * whole frame rather than three slots.  Same abort contract as vm_array_op. */
/* Values[assoc] -> owned packed vector (B1).  The source is the spec's constant
 * association or, for an argument bag, the borrowed handle in R[c->a].  Declines
 * (returns false -> interpreter) on a non-numeric / ragged value list. */
static bool vm_assoc_values(const Instr* c, Slot* R) {
    const AssocSpec* sp = (const AssocSpec*)c->imm.p;
    Expr* assoc = sp->assoc ? sp->assoc : R[c->a].arr;
    if (!assoc) return false;
    Expr* vals = assoc_values_list(assoc);        /* owned List of the values */
    if (!vals) return false;
    Expr* nd = ndarray_from_nested_list(vals, assoc_elem_ndt((CompileType)(c->flags & 0xFFu)));
    expr_free(vals);
    if (!nd) return false;                         /* non-numeric -> decline    */
    if ((c->flags & 0x100u) && !sp->assoc) { expr_free(R[c->a].arr); R[c->a].arr = NULL; }  /* free produced src */
    expr_free(R[c->dst].arr);                      /* release any stale handle  */
    R[c->dst].arr = nd;
    return true;
}

/* KeyDrop/KeyTake (B3) -> an OWNED association (array bank).  Native: builds the
 * filtered association directly (assoc_key_select), no evaluator, no call node.
 * The result carries Part A's index so a downstream read stays O(1).  flags bit0
 * = take; bit1 = free the source temp (a produced association consumed here). */
static bool vm_assoc_keysel(const Instr* c, Slot* R) {
    const AssocSpec* sp = (const AssocSpec*)c->imm.p;
    Expr* src = sp->assoc ? sp->assoc : R[c->a].arr;
    if (!src) return false;
    Expr* r = assoc_key_select(src, sp->key, (c->flags & 1u) != 0);
    if (!r || !is_association(r)) { expr_free(r); return false; }
    assoc_prebuild_index(r);
    if ((c->flags & 2u) && !sp->assoc) { expr_free(R[c->a].arr); R[c->a].arr = NULL; }
    expr_free(R[c->dst].arr);
    R[c->dst].arr = r;
    return true;
}

static bool vm_call(const CompiledProgram* cp, const Slot* argv, unsigned nargs, Slot* dst);

/* B4: Map[f, assoc] -> a new OWNED association, same keys, each value replaced by
 * the compiled callee's result.  Select[assoc, pred] keeps the entries whose value
 * the Bool callee accepts.  The callee runs per entry via vm_call — machine value
 * in, machine result out, NO evaluator.  flags = in_valtype | out_type<<4 |
 * free-source<<8. */
static bool vm_assoc_higher(const Instr* c, Slot* R, bool select) {
    const AssocCalleeSpec* sp = (const AssocCalleeSpec*)c->imm.p;
    Expr* src = sp->assoc ? sp->assoc : R[c->a].arr;
    if (!src || !is_association(src)) return false;
    CompileType in_t  = (CompileType)(c->flags & 0xFu);
    CompileType out_t = (CompileType)((c->flags >> 4) & 0xFu);
    size_t n = src->data.function.arg_count;
    Expr** out = malloc((n ? n : 1) * sizeof(Expr*));
    if (!out) return false;
    size_t nout = 0;
    bool ok = true;
    for (size_t i = 0; i < n; i++) {
        Expr* entry = src->data.function.args[i];
        if (!is_rule2(entry)) { ok = false; break; }
        Slot arg;
        if (!assoc_value_to_slot(entry->data.function.args[1], in_t, &arg)) { ok = false; break; }
        Slot res;
        if (!vm_call(sp->callee, &arg, 1, &res)) { ok = false; break; }
        if (select) {
            if (res.i) out[nout++] = expr_copy(entry);            /* keep entry     */
        } else {
            Expr* nv = assoc_slot_to_value(res, out_t);
            if (!nv) { ok = false; break; }
            out[nout++] = assoc_entry_with_value(entry, nv);      /* same key, new v */
        }
    }
    if (!ok) { for (size_t j = 0; j < nout; j++) expr_free(out[j]); free(out); return false; }
    Expr* result = expr_new_function(expr_new_symbol(SYM_Association), out, nout);
    free(out);
    if (!result) return false;
    assoc_prebuild_index(result);
    if ((c->flags & 0x100u) && !sp->assoc) { expr_free(R[c->a].arr); R[c->a].arr = NULL; }
    expr_free(R[c->dst].arr);
    R[c->dst].arr = result;
    return true;
}

/* Append[assoc, key -> value] (B5) -> a new OWNED association with key set
 * (replaced in place, else appended).  Native (assoc_set_key), no evaluator.
 * R[a] = source, R[b] = machine value; flags = value_type | free-source<<8. */
static bool vm_assoc_set(const Instr* c, Slot* R) {
    const AssocSpec* sp = (const AssocSpec*)c->imm.p;
    Expr* src = sp->assoc ? sp->assoc : R[c->a].arr;
    if (!src || !is_association(src)) return false;
    Expr* nv = assoc_slot_to_value(R[c->b], (CompileType)(c->flags & 0xFu));
    if (!nv) return false;
    Expr* r = assoc_set_key(src, sp->key, nv);      /* adopts nv */
    if (!r) { expr_free(nv); return false; }
    assoc_prebuild_index(r);
    if ((c->flags & 0x100u) && !sp->assoc) { expr_free(R[c->a].arr); R[c->a].arr = NULL; }
    expr_free(R[c->dst].arr);
    R[c->dst].arr = r;
    return true;
}

/* Counts[machine array] (B3) -> an OWNED association of element->count.  Native
 * (assoc_counts_ndarray via Tally, no evaluator).  flags bit0 = free the source
 * array temp (a produced array consumed here). */
static bool vm_assoc_counts(const Instr* c, Slot* R) {
    Expr* arr = R[c->a].arr;
    if (!arr) return false;
    Expr* r = assoc_counts_ndarray(arr);
    if (!r || !is_association(r)) { expr_free(r); return false; }
    assoc_prebuild_index(r);
    if (c->flags & 1u) { expr_free(R[c->a].arr); R[c->a].arr = NULL; }
    expr_free(R[c->dst].arr);
    R[c->dst].arr = r;
    return true;
}

static bool vm_range_array_op(const Instr* c, Slot* R) {
    switch (c->op) {
        case OP_A_NEW: {                  /* ConstantArray[0, {d1, ..., dr}] */
            int rank = (int)c->flags;
            if (rank < 1 || rank > NDARRAY_MAX_RANK) return false;
            int64_t dims[NDARRAY_MAX_RANK];
            size_t n = 1;
            for (int i = 0; i < rank; i++) {
                long long d = R[c->a + (unsigned)i].i;
                if (d < 0) return false;          /* Table[..., {n}] with n < 0 */
                dims[i] = (int64_t)d;
                n *= (size_t)d;
            }
            NDType dt = ct_elem_ndt((unsigned)c->imm.i);
            void* buf = calloc(n ? n : 1, ndt_elem_size(dt));
            if (!buf) return false;
            Expr* nw = expr_new_ndarray_raw(rank, dims, buf, dt);
            if (!nw) { free(buf); return false; }
            expr_free(R[c->dst].arr);             /* register reused from a prior call */
            R[c->dst].arr = nw;
            return true;
        }

        case OP_A_PART: {                 /* Span / All / list / partial indexing */
            const PartSpec* ps = (const PartSpec*)c->imm.p;
            Expr* idx[NDARRAY_MAX_RANK];
            const Expr* src = R[c->b].arr;
            if (!src || src->type != EXPR_NDARRAY) return false;
            if (!ps_build_indices(ps, R, idx)) return false;
            bool degrade = false;
            Expr* r = ndarray_part(src, idx, (size_t)ps->n, &degrade);
            ps_free_indices(ps, idx);
            /* A spec ndarray_part cannot do natively (degrade) or an out-of-range
             * subscript (r == NULL) both abort: the interpreter re-runs the body
             * and produces whatever Part[] properly produces, including a
             * diagnostic.  The compiled path never invents an answer. */
            if (!r) return false;
            if (r->type != EXPR_NDARRAY) { expr_free(r); return false; }
            expr_free(R[c->dst].arr);
            R[c->dst].arr = r;
            return true;
        }

        case OP_A_NDFN: {                 /* Reverse / Sort / Flatten / Take / ... */
            const NdFnSpec* fn = (const NdFnSpec*)c->imm.p;
            Expr* src = R[c->b].arr;
            if (!src || src->type != EXPR_NDARRAY) return false;
            size_t nx = (size_t)c->flags;
            /* The entry point takes the whole CALL, so rebuild it.  expr_copy is
             * a refcount bump (src/expr.c), not a buffer copy, so this costs a
             * node — and the array is never mutated: these paths all allocate a
             * fresh result. */
            Expr* args[1 + NDARRAY_MAX_RANK];
            args[0] = expr_copy(src);
            for (size_t i = 0; i < nx; i++)
                args[1 + i] = expr_new_integer((int64_t)R[c->a + (unsigned)i].i);
            Expr* call = expr_new_function(expr_new_symbol(fn->head), args, 1 + nx);
            if (!call) { expr_free(args[0]); return false; }
            Expr* r = fn->fn(call);
            expr_free(call);
            /* A case the fast path does not handle comes back as a nested List
             * (ndarray_delist_and_reeval), which is precisely the signal to
             * decline: the interpreter then answers, exactly as it would have. */
            if (!r) return false;
            if (r->type != EXPR_NDARRAY) { expr_free(r); return false; }
            expr_free(R[c->dst].arr);
            R[c->dst].arr = r;
            return true;
        }

        case OP_A_PARTSET: {              /* u[[spec...]] = rhs, in place */
            const PartSpec* ps = (const PartSpec*)c->imm.p;
            Expr* idx[NDARRAY_MAX_RANK];
            Expr* tgt = R[c->dst].arr;
            if (!tgt || tgt->type != EXPR_NDARRAY) return false;
            if (!ps_build_indices(ps, R, idx)) return false;
            Expr* rhs = (ps->rhs_kind == AK_ARR)
                      ? R[c->b].arr
                      : vm_box_scalar(&R[c->b], (unsigned)ps->rhs_kind);
            bool ok = rhs && ndarray_part_set(tgt, idx, (size_t)ps->n, rhs);
            if (ps->rhs_kind != AK_ARR) expr_free(rhs);
            ps_free_indices(ps, idx);
            return ok;
        }

        default: return false;
    }
}

static bool vm_array_op(const Instr* c, Slot* d, Slot* a, Slot* b) {
    const unsigned f = c->flags, ka = AF_A(f), kb = AF_B(f);
    Expr* r = NULL;

    switch (c->op) {
        case OP_ARR_FREE:
            expr_free(d->arr); d->arr = NULL;
            return true;

        /* ---- ownership (M3c) ------------------------------------------ */
        case OP_A_COPY: {                 /* a local initialised from a borrowed array */
            Expr* nw = nd_own_copy(a->arr, AF_R(f));
            if (!nw) return false;
            if (f & AF_FREE_A) { expr_free(a->arr); a->arr = NULL; }
            expr_free(d->arr);
            d->arr = nw;
            return true;
        }

        case OP_A_XFER:                   /* move the handle, do not duplicate it */
            if (d == a) return true;
            expr_free(d->arr);
            d->arr = a->arr;
            a->arr = NULL;                /* or teardown would free it twice */
            return true;

        /* ---- fused-loop setup (M3b) -----------------------------------
         * Out of line because each runs ONCE per call, outside the element
         * loop; only A_LOAD/A_STORE are inline in the dispatch loop. */
        case OP_A_SIZE: {                 /* total elements, any rank */
            const Expr* x = a->arr;
            if (!x || x->type != EXPR_NDARRAY) return false;
            d->i = (long long)ndarray_size(x);
            return true;
        }

        case OP_A_SHAPECHK: {             /* all leaves of a fused loop agree */
            const Expr* x = a->arr;
            const Expr* y = b->arr;
            if (!x || !y || x->type != EXPR_NDARRAY || y->type != EXPR_NDARRAY) return false;
            if (x->data.ndarray.rank != y->data.ndarray.rank) return false;
            for (int i = 0; i < x->data.ndarray.rank; i++)
                if (x->data.ndarray.dims[i] != y->data.ndarray.dims[i]) return false;
            return true;
        }

        case OP_A_NEWLIKE: {              /* result buffer, shape of the operand */
            const Expr* x = a->arr;
            if (!x || x->type != EXPR_NDARRAY) return false;
            NDType dt = ct_elem_ndt(AF_R(f));
            size_t nelem = ndarray_size(x);
            size_t esz = ndt_elem_size(dt);
            void* buf = calloc(nelem ? nelem : 1, esz);
            if (!buf) return false;
            Expr* nw = expr_new_ndarray_like(x, x->data.ndarray.rank, x->data.ndarray.dims, buf, dt);
            if (!nw) { free(buf); return false; }
            expr_free(d->arr);            /* reused register from a prior call */
            d->arr = nw;
            return true;
        }

        case OP_V_LEN: {
            const Expr* x = a->arr;
            if (!x || x->type != EXPR_NDARRAY) return false;
            long long len = (long long)x->data.ndarray.dims[0];
            if (f & AF_FREE_A) { expr_free(a->arr); a->arr = NULL; }
            d->i = len;
            return true;
        }

        case OP_V_TOTAL: {
            /* An INTEGER array sums here rather than through ndred_total_all,
             * which accumulates in a double and is therefore exact only to 2^53
             * — `Total[{9007199254740993, 1}]` came back one short.  Summing in
             * int64 with the same overflow rule as the scalar opcodes keeps the
             * answer identical to the interpreter's, or hands the call back. */
            if (a->arr && a->arr->type == EXPR_NDARRAY
                && a->arr->data.ndarray.dtype == NDT_INT64) {
                const NDArrayData* A_ = &a->arr->data.ndarray;
                const int64_t* p = (const int64_t*)A_->data;
                size_t n = ndarray_size(a->arr);
                long long acc = 0;
                bool ovf = false;
                for (size_t k = 0; k < n && !ovf; k++) ovf = ci_add(acc, (long long)p[k], &acc);
                if (f & AF_FREE_A) { expr_free(a->arr); a->arr = NULL; }
                if (ovf) return false;
                d->i = acc;
                return true;
            }
            Expr* s = ndred_total_all(a->arr);       /* borrows; same rounding as Total[] */
            if (f & AF_FREE_A) { expr_free(a->arr); a->arr = NULL; }
            bool ok = s && vm_write_scalar(s, AF_R(f), d);
            expr_free(s);
            return ok;
        }

        case OP_V_NDRED: {                /* Mean / Median / Variance / Max / ... */
            const NdRedSpec* rs = (const NdRedSpec*)c->imm.p;
            if (!a->arr || a->arr->type != EXPR_NDARRAY) return false;
            /* The entry point takes the whole CALL; expr_copy is a refcount
             * bump (src/expr.c), so rebuilding it costs two nodes and never
             * the buffer.  Same shape as OP_A_NDFN above. */
            Expr* arg  = expr_copy(a->arr);
            Expr* call = expr_new_function(expr_new_symbol(rs->head), &arg, 1);
            if (!call) { expr_free(arg); return false; }
            Expr* s = rs->fn(call);
            expr_free(call);
            if (f & AF_FREE_A) { expr_free(a->arr); a->arr = NULL; }
            /* A form the reduction declines comes back as a materialised List
             * (ndarray_delist_and_reeval) rather than a number, and
             * vm_write_scalar refuses it -- which aborts the program and lets
             * the interpreter answer, exactly as it would have. */
            bool ok = s && vm_write_scalar(s, AF_R(f), d);
            expr_free(s);
            return ok;
        }

        case OP_V_NDREDN: {               /* RankedMin[v, n] / RankedMax[v, n] */
            const NdRedSpec* rs = (const NdRedSpec*)c->imm.p;
            if (!a->arr || a->arr->type != EXPR_NDARRAY) return false;
            /* V_NDRED's rebuild-and-delegate, plus one trailing integer read from
             * the b register (A_NDFN's trick).  A runtime-out-of-range n comes
             * back as a materialised List, which vm_write_scalar refuses -> the
             * interpreter answers, exactly as it would have. */
            Expr* args[2] = { expr_copy(a->arr), expr_new_integer((int64_t)b->i) };
            Expr* call = expr_new_function(expr_new_symbol(rs->head), args, 2);
            if (!call) { expr_free(args[0]); expr_free(args[1]); return false; }
            Expr* s = rs->fn(call);
            expr_free(call);
            if (f & AF_FREE_A) { expr_free(a->arr); a->arr = NULL; }
            bool ok = s && vm_write_scalar(s, AF_R(f), d);
            expr_free(s);
            return ok;
        }

        case OP_A_NDFN2: {                /* Dot (matrix) / LinearSolve / Cross / Join / ... */
            const NdFn2Spec* fn = (const NdFn2Spec*)c->imm.p;
            if (!a->arr || a->arr->type != EXPR_NDARRAY
                || !b->arr || b->arr->type != EXPR_NDARRAY) return false;
            /* Rebuild the whole call and delegate, exactly as A_NDFN — expr_copy
             * is a refcount bump, so this costs two nodes and never a buffer. */
            Expr* args2[2] = { expr_copy(a->arr), expr_copy(b->arr) };
            Expr* call = expr_new_function(expr_new_symbol(fn->head), args2, 2);
            if (!call) { expr_free(args2[0]); expr_free(args2[1]); return false; }
            r = fn->fn(call);
            expr_free(call);
            /* Positive result-dtype guard: the promised element type (AF_R) must
             * match the buffer the delegate returned, or a downstream fused loop
             * would read it at the wrong element width.  A mismatch — or a
             * declined non-NDArray (a scalar vector.vector Dot, a nested List)
             * — falls through the shared tail below to the interpreter. */
            if (r && r->type == EXPR_NDARRAY
                && r->data.ndarray.dtype != ct_elem_ndt(AF_R(f))) {
                expr_free(r); r = NULL;
            }
            break;   /* -> shared array tail: frees both operands, stores or declines */
        }

        case OP_V_NDFN2: {                /* Dot's vector.vector inner product -> scalar */
            const NdFn2Spec* fn = (const NdFn2Spec*)c->imm.p;
            if (!a->arr || a->arr->type != EXPR_NDARRAY
                || !b->arr || b->arr->type != EXPR_NDARRAY) return false;
            Expr* args2[2] = { expr_copy(a->arr), expr_copy(b->arr) };
            Expr* call = expr_new_function(expr_new_symbol(fn->head), args2, 2);
            if (!call) { expr_free(args2[0]); expr_free(args2[1]); return false; }
            Expr* s = fn->fn(call);
            expr_free(call);
            /* BOTH operands are freed here — V_NDRED, its unary template, frees
             * only one; a binary op that dropped the second would leak it. */
            if (f & AF_FREE_A) { expr_free(a->arr); a->arr = NULL; }
            if (f & AF_FREE_B) { expr_free(b->arr); b->arr = NULL; }
            bool ok = s && vm_write_scalar(s, AF_R(f), d);
            expr_free(s);
            return ok;
        }

        case OP_V_EW: {
            Expr* boxa = (ka == AK_ARR) ? NULL : vm_box_scalar(a, ka);
            Expr* boxb = (kb == AK_ARR) ? NULL : vm_box_scalar(b, kb);
            Expr* ops[2];
            ops[0] = boxa ? boxa : a->arr;
            ops[1] = boxb ? boxb : b->arr;
            if (ops[0] && ops[1]) r = ndarray_elementwise(ops, 2, c->imm.i != 0);
            expr_free(boxa); expr_free(boxb);
            break;
        }

        case OP_V_POW: {
            double re, im;
            if (ka == AK_ARR && kb == AK_ARR)  r = ndarray_elementwise_power(a->arr, b->arr);
            else if (ka == AK_ARR) { vm_scalar_pair(b, kb, &re, &im); r = ndarray_scalar_power(a->arr, re, im); }
            else                   { vm_scalar_pair(a, ka, &re, &im); r = ndarray_base_scalar_power(re, im, b->arr); }
            break;
        }

        case OP_V_KERN:
            r = ndarray_map_unary(a->arr, (const NDUnaryKernel*)c->imm.p);
            break;

        case OP_V_KERN2: {
            Expr* boxa = (ka == AK_ARR) ? NULL : vm_box_scalar(a, ka);
            Expr* boxb = (kb == AK_ARR) ? NULL : vm_box_scalar(b, kb);
            const Expr* p0 = boxa ? boxa : a->arr;
            const Expr* p1 = boxb ? boxb : b->arr;
            if (p0 && p1) r = ndarray_map_binary(p0, p1, (const NDBinaryKernel*)c->imm.p);
            expr_free(boxa); expr_free(boxb);
            break;
        }

        default: return false;
    }

    /* shared tail for the array-PRODUCING ops */
    if (f & AF_FREE_A) { expr_free(a->arr); a->arr = NULL; }
    if (f & AF_FREE_B) { expr_free(b->arr); b->arr = NULL; }
    if (!r) return false;
    if (r->type != EXPR_NDARRAY
        || (AF_R(f) != (unsigned)CT_COMPLEX && ndt_is_complex(r->data.ndarray.dtype))) {
        expr_free(r); return false;
    }
    d->arr = r;
    return true;
}

/* Every opcode, for the computed-goto jump table (must cover the whole enum). */
/* OPLIST now lives in compile_internal.h: one list drives the opcode enum, the
 * VM jump table below, and the optimiser's instruction-property table. */

/* The bytecode interpreter.  A threaded (computed-goto) dispatch is used on
 * GCC/Clang — each opcode ends by jumping straight to the next, which the branch
 * predictor handles far better than a single switch; a portable switch is the
 * fallback.  Programs always end in OP_RET, so no per-op bounds check is needed. */
#if defined(__GNUC__) && !defined(VM_NO_THREADED)
#define VM_THREADED 1
#else
#define VM_THREADED 0
#endif

static bool vm_call(const CompiledProgram* cp, const Slot* argv, unsigned nargs, Slot* dst);
static void vm_run(const Instr* code, size_t n, Slot* R, bool* failed);

/* ------------------------------------------------------------------ *
 *  Parallel strip loop (OP_APAR)                                      *
 * ------------------------------------------------------------------ *
 * A fused MAP is embarrassingly parallel: element i of the output depends only
 * on element i of the inputs, so a worker can be handed a sub-range of the flat
 * index space and the result is BIT-IDENTICAL to the serial pass — same
 * operations, same order, on the same values. Nothing here changes an answer;
 * it only changes which core computes it.
 *
 * Each worker gets its OWN frame, seeded by copying the parent's registers.
 * That copy is what makes this safe with no locking anywhere:
 *   - scalars and loop bounds are inherited by value;
 *   - array registers are inherited as raw pointers, which is correct because
 *     the workers only READ the inputs and write DISJOINT output elements;
 *   - tile registers are re-pointed at the worker's own tile storage, so no two
 *     workers ever touch the same VBLOCK scratch buffer.
 * No worker allocates or frees an array (the output buffer was allocated by the
 * parent before the loop), so no ownership crosses a thread boundary.
 *
 * Returns false to DECLINE — too small to be worth threading, out of memory, or
 * a worker failed. Declining is always safe: OP_APAR then falls through into the
 * serial loop, which is still sitting immediately after it and which recomputes
 * the same map from the same starting index. */
#define VM_STACK_SLOTS 512

typedef struct { const ParLoop* pl; const Slot* parent; } ParCtx;

static bool par_chunk(void* vctx, size_t lo, size_t hi) {
    const ParCtx* x = (const ParCtx*)vctx;
    const ParLoop* pl = x->pl;

    Slot stackframe[VM_STACK_SLOTS];
    Slot* heap = NULL;
    Slot* R = stackframe;
    if (pl->frame_slots > VM_STACK_SLOTS) {
        heap = malloc(pl->frame_slots * sizeof(Slot));
        if (!heap) return false;
        R = heap;
    }
    memcpy(R, x->parent, (size_t)pl->nreg * sizeof(Slot));
    Slot* tiles = R + pl->nreg;
    for (int k = 0; k < pl->ntiles; k++)
        R[pl->tile_base + k].p = tiles + (size_t)k * VBLOCK;

    R[pl->ri].i = (long long)lo;
    R[pl->rn].i = (long long)hi;


    bool failed = false;
    vm_run(pl->code, pl->n, R, &failed);
    free(heap);
    return !failed;
}

static bool vm_par_run(const ParLoop* pl, Slot* R) {
#ifdef MATHILDA_THREADS
    long long n = R[pl->rn].i;
    /* nd_parallel_for owns the "is this worth threading" decision and the
     * chunking, so the compiled path fans out on exactly the same terms as the
     * rest of the ND layer rather than inventing a second policy. It runs the
     * range serially below its threshold, which for us would be pure overhead —
     * hence the explicit check first. */
    if (n < (long long)NDARRAY_THREAD_THRESHOLD) return false;

    ParCtx ctx; ctx.pl = pl; ctx.parent = R;
    return nd_parallel_for((size_t)n, par_chunk, &ctx);
#else
    (void)pl; (void)R;
    return false;
#endif
}

static void vm_run(const Instr* code, size_t n, Slot* R, bool* failed) {
    *failed = false;
    if (n == 0) return;
    size_t pc = 0;
    const Instr* c = &code[pc];
    /* Live elements of the current tile, set by VSETLEN.  Only the load, store
     * and reduction read it — every arithmetic tile op covers the full VBLOCK —
     * and all three are pinned inside the strip loop because they read the loop
     * index, so this cannot be stranded by a hoist. */
    int vlen = VBLOCK;
    /* Operands are addressed lazily: an opcode pays only for the registers it
     * actually reads.  Computing all three up front in NEXT() cost three
     * shift-and-adds per instruction where the hottest ops (ADD_R, MUL_R) use
     * two and the unary ops use one. */
    #define RD (R[c->dst])
    #define RA (R[c->a])
    #define RB (R[c->b])
#if VM_THREADED
    #define X(name, kind) [OP_##name] = &&L_##name,
    static const void* const tbl[] = { OPLIST };
    #undef X
    #define OP(name) L_##name
    #define NEXT() do { c = &code[++pc]; goto *tbl[c->op]; } while (0)
    #define JUMP() do { c = &code[pc]; goto *tbl[c->op]; } while (0)
    goto *tbl[c->op];
#else
    #define OP(name) case OP_##name
    #define NEXT() break
    #define JUMP() continue
    while (pc < n) {
        c = &code[pc];
        switch (c->op) {
#endif
    #define ARROP() do { if (!vm_array_op(c, &RD, &RA, &RB)) goto vm_fail; } while (0); NEXT()
    /* An integer opcode that can overflow.  The helper returns true on overflow,
     * and the call is abandoned exactly as OP_FAIL abandons it, so the caller
     * falls back to the interpreter and gets the exact bigint answer.
     *
     * Two levels of opting out, and they measure different things:
     *
     *   IF_NOCHK in the instruction (COMPILE_WRAP_INT / Compile[]'s
     *   "CatchMachineIntegerOverflow" -> False) keeps the wrapped value instead
     *   of falling back.  Because `&&` short-circuits, the flag is read only
     *   once an overflow has actually been detected, so a program that does not
     *   overflow executes the identical instruction path either way — the option
     *   costs nothing.
     *
     *   `-DVM_NO_INT_CHECK` at BUILD time removes the detection itself.  That is
     *   the one that measures what the feature really costs, because it is the
     *   detection — not the never-taken branch — that is the price.  (Named for
     *   the VM, like VM_NO_THREADED, and deliberately NOT after the
     *   COMPILE_WRAP_INT flag: compile.h defines that as a bit value, so an
     *   `#ifdef` on it is true in every build and silently disables the
     *   checks — which is exactly the bug this comment now prevents.) */
#ifdef VM_NO_INT_CHECK
    #define IOP(chk) do { (void)(chk); } while (0); NEXT()
#else
    #define IOP(chk) do { if ((chk) && !(c->flags & IF_NOCHK)) goto vm_fail; } while (0); NEXT()
#endif
            OP(JMP): pc = c->b; JUMP();
            OP(JZ):  pc = RA.i ? pc + 1 : c->b; JUMP();   /* branch if false */
            /* Checked: this is Increment/Decrement's opcode as well as a loop
             * step (compile.c's Increment lowering), and `x = 2^63-1; x++` is a
             * perfectly ordinary thing to write. */
            OP(INC_I): IOP(ci_add(RD.i, c->imm.i, &RD.i));
            /* Increment, test and branch in one.  A counted loop otherwise spends
             * four instructions per iteration on control alone (INC, LT, JZ,
             * JMP), which is most of the body when the body is one element of a
             * fused array pass.
             *
             * Deliberately NOT overflow-checked, unlike INC_I.  The index only
             * overflows if the limit sits within one step of INT64_MAX, and
             * reaching there from the initial value takes on the order of 10^18
             * iterations — no terminating run gets close.  The check would
             * otherwise land in the innermost loop of every fused array pass. */
            OP(LOOP): { RD.i += c->imm.i; if (RD.i < RA.i) { pc = c->b; JUMP(); } } NEXT();
            /* Run the strip loop that follows across threads and skip past it;
             * on a decline, fall through and run it right here, serially. */
            OP(APAR): if (vm_par_run((const ParLoop*)c->imm.p, R)) { pc = c->b; JUMP(); }
                      NEXT();
            OP(CONST): RD = c->imm; NEXT();
            OP(MOVE):  RD = RA; NEXT();
            OP(I2R): RD.r = (double)RA.i; NEXT();
            OP(I2C): RD.z = (double)RA.i; NEXT();
            OP(R2C): RD.z = RA.r; NEXT();
            OP(ADD_I): IOP(ci_add(RA.i, RB.i, &RD.i));
            OP(ADD_R): RD.r = RA.r + RB.r; NEXT();
            OP(ADD_C): RD.z = RA.z + RB.z; NEXT();
            /* Immediate forms (K_BINK).  One register read instead of two, and
             * the CONST that would have materialised the operand is gone.  Each
             * is a SINGLE arithmetic operation, so no floating-point contraction
             * is possible and the result is bit-identical to the register form
             * the optimiser rewrote. */
            OP(ADD_RK): RD.r = RA.r + c->imm.r; NEXT();
            OP(SUB_RK): RD.r = RA.r - c->imm.r; NEXT();
            OP(SUB_KR): RD.r = c->imm.r - RA.r; NEXT();
            OP(MUL_RK): RD.r = RA.r * c->imm.r; NEXT();
            OP(DIV_RK): RD.r = RA.r / c->imm.r; NEXT();
            OP(DIV_KR): RD.r = c->imm.r / RA.r; NEXT();
            OP(ADD_IK): IOP(ci_add(RA.i, c->imm.i, &RD.i));
            OP(SUB_IK): IOP(ci_sub(RA.i, c->imm.i, &RD.i));
            OP(SUB_KI): IOP(ci_sub(c->imm.i, RA.i, &RD.i));
            OP(MUL_IK): IOP(ci_mul(RA.i, c->imm.i, &RD.i));
            OP(LT_RK): RD.i = RA.r <  c->imm.r; NEXT();
            OP(LE_RK): RD.i = RA.r <= c->imm.r; NEXT();
            OP(GT_RK): RD.i = RA.r >  c->imm.r; NEXT();
            OP(GE_RK): RD.i = RA.r >= c->imm.r; NEXT();
            OP(LT_IK): RD.i = RA.i <  c->imm.i; NEXT();
            OP(LE_IK): RD.i = RA.i <= c->imm.i; NEXT();
            OP(GT_IK): RD.i = RA.i >  c->imm.i; NEXT();
            OP(GE_IK): RD.i = RA.i >= c->imm.i; NEXT();
            OP(SUB_I): IOP(ci_sub(RA.i, RB.i, &RD.i));
            OP(SUB_R): RD.r = RA.r - RB.r; NEXT();
            OP(SUB_C): RD.z = RA.z - RB.z; NEXT();
            OP(MUL_I): IOP(ci_mul(RA.i, RB.i, &RD.i));
            OP(MUL_R): RD.r = RA.r * RB.r; NEXT();
            OP(MUL_C): RD.z = RA.z * RB.z; NEXT();
            OP(DIV_R): RD.r = RA.r / RB.r; NEXT();
            OP(DIV_C): RD.z = RA.z / RB.z; NEXT();
            /* Both integer divisions guard the two inputs the hardware traps on
             * rather than merely answering wrongly: a zero divisor, and
             * INT64_MIN / -1 whose quotient is not representable.  `Mod[5, 0]`
             * is left unevaluated by the interpreter, and the compiled path used
             * to take the whole process down with SIGFPE. */
            OP(MOD_I): { long long m = RB.i, x = RA.i;
                         if (m == 0 || (m == -1 && x == LLONG_MIN)) goto vm_fail;
                         long long q = x % m; if (q != 0 && ((q < 0) != (m < 0))) q += m;
                         RD.i = q; } NEXT();
            OP(QUOT_I): { long long m = RB.i, x = RA.i;
                          if (m == 0 || (m == -1 && x == LLONG_MIN)) goto vm_fail;
                          long long q = x / m; if ((x % m != 0) && ((x < 0) != (m < 0))) q -= 1;
                          RD.i = q; } NEXT();
            OP(NEG_I): IOP(ci_neg(RA.i, &RD.i));
            OP(NEG_R): RD.r = -RA.r; NEXT();
            OP(NEG_C): RD.z = -RA.z; NEXT();
            OP(INV_R): RD.r = 1.0 / RA.r; NEXT();
            OP(INV_C): RD.z = 1.0 / RA.z; NEXT();
            OP(POWI_I): IOP(ci_powi(RA.i, c->imm.i, &RD.i));
            OP(POWI_R): RD.r = ipow_r(RA.r, c->imm.i); NEXT();
            OP(POWI_C): RD.z = ipow_c(RA.z, c->imm.i); NEXT();
            OP(POW_R): RD.r = pow(RA.r, RB.r); NEXT();
            OP(POW_C): RD.z = cpow(RA.z, RB.z); NEXT();
            OP(SQRT_R): RD.r = sqrt(RA.r); NEXT();
            OP(SQRT_C): RD.z = csqrt(RA.z); NEXT();
            OP(EXP_R): RD.r = exp(RA.r); NEXT();
            OP(EXP_C): RD.z = cexp(RA.z); NEXT();
            OP(LOG_R): RD.r = log(RA.r); NEXT();
            OP(LOG_C): RD.z = clog(RA.z); NEXT();
            OP(SIN_R): RD.r = sin(RA.r); NEXT();   OP(SIN_C): RD.z = csin(RA.z); NEXT();
            OP(COS_R): RD.r = cos(RA.r); NEXT();   OP(COS_C): RD.z = ccos(RA.z); NEXT();
            OP(TAN_R): RD.r = tan(RA.r); NEXT();   OP(TAN_C): RD.z = ctan(RA.z); NEXT();
            OP(SINH_R): RD.r = sinh(RA.r); NEXT(); OP(SINH_C): RD.z = csinh(RA.z); NEXT();
            OP(COSH_R): RD.r = cosh(RA.r); NEXT(); OP(COSH_C): RD.z = ccosh(RA.z); NEXT();
            OP(TANH_R): RD.r = tanh(RA.r); NEXT(); OP(TANH_C): RD.z = ctanh(RA.z); NEXT();
            OP(ASIN_R): RD.r = asin(RA.r); NEXT(); OP(ASIN_C): RD.z = casin(RA.z); NEXT();
            OP(ACOS_R): RD.r = acos(RA.r); NEXT(); OP(ACOS_C): RD.z = cacos(RA.z); NEXT();
            OP(ATAN_R): RD.r = atan(RA.r); NEXT(); OP(ATAN_C): RD.z = catan(RA.z); NEXT();
            OP(ABS_I): IOP(ci_abs(RA.i, &RD.i));
            OP(ABS_R): RD.r = fabs(RA.r); NEXT();
            OP(ABS_C): RD.r = cabs(RA.z); NEXT();
            OP(SIGN_I): RD.i = (RA.i > 0) - (RA.i < 0); NEXT();
            /* Integer-closed heads.  Same IOP() shape as the arithmetic
             * opcodes: the helper reports overflow OR an out-of-domain
             * argument, and either way the interpreter takes the call. */
            OP(POW_II):  IOP(int_pow(RA.i, RB.i, &RD.i));
            OP(FACT_I):  IOP(int_factorial(RA.i, &RD.i));
            OP(GAMMA_I): IOP(int_gamma(RA.i, &RD.i));
            OP(BINOM_I): IOP(int_binomial(RA.i, RB.i, &RD.i));
            OP(POCH_I):  IOP(int_pochhammer(RA.i, RB.i, &RD.i));
            OP(FIB_I):   IOP(int_fib(RA.i, &RD.i));
            OP(LUCAS_I): IOP(int_lucas(RA.i, &RD.i));
            /* Arg[n] is 0 for n >= 0; below that it is the symbol Pi. */
            OP(ARG_I):   { if (RA.i < 0) goto vm_fail; RD.i = 0; } NEXT();
            OP(GCD_I):   IOP(int_gcd(RA.i, RB.i, &RD.i));
            OP(LCM_I):   IOP(int_lcm(RA.i, RB.i, &RD.i));
            OP(ILEN_I):  IOP(int_ilen(RA.i, RB.i, &RD.i));
            OP(IEXP_I):  IOP(int_iexp(RA.i, RB.i, &RD.i));
            /* Ternary, so the operands are a RUN of registers starting at `a`
             * (the K_NARY shape KERNN uses), not the two operand fields. */
            OP(POWMOD_I): IOP(int_powmod(R[c->a].i, R[c->a + 1].i, R[c->a + 2].i, &RD.i));
            /* Writes the INTEGER slot: Sign of a real is an Integer in the
             * interpreter.  (The tile form VSIGN_R stays real — it feeds packed
             * float arrays, where the ND kernel's result dtype is real.) */
            OP(SIGN_R): RD.i = (RA.r > 0) - (RA.r < 0); NEXT();
            OP(FLOOR_R): RD.i = (long long)floor(RA.r); NEXT();
            /* IntegerPart truncates TOWARD ZERO, which is not Floor for a
             * negative argument: IntegerPart[-1.5] is -1, Floor[-1.5] is -2. */
            OP(TRUNC_R): RD.i = (long long)trunc(RA.r); NEXT();
            OP(CEIL_R):  RD.i = (long long)ceil(RA.r); NEXT();
            OP(ROUND_R): RD.i = (long long)llround(RA.r); NEXT();
            OP(RE_C): RD.r = creal(RA.z); NEXT();
            OP(IM_C): RD.r = cimag(RA.z); NEXT();
            OP(ARG_C): RD.r = carg(RA.z); NEXT();
            OP(CONJ_C): RD.z = conj(RA.z); NEXT();
            OP(ATAN2_R): RD.r = atan2(RB.r, RA.r); NEXT();   /* ArcTan[x,y]=atan2(y,x); a=x,b=y */
            OP(MAX_I): RD.i = RA.i > RB.i ? RA.i : RB.i; NEXT();
            OP(MAX_R): RD.r = RA.r > RB.r ? RA.r : RB.r; NEXT();
            OP(MIN_I): RD.i = RA.i < RB.i ? RA.i : RB.i; NEXT();
            OP(MIN_R): RD.r = RA.r < RB.r ? RA.r : RB.r; NEXT();
            OP(ERF_R): RD.r = erf(RA.r); NEXT();
            OP(ERFC_R): RD.r = erfc(RA.r); NEXT();
            OP(KERN_RR): { double o; RD.r = ((kfn_r)c->imm.p)(RA.r, &o) ? o : NAN; } NEXT();
            OP(KERN_R2R): { double orr, oi; RD.r = ((kfn_c)c->imm.p)(RA.r, 0.0, &orr, &oi) ? orr : NAN; } NEXT();
            OP(KERN_RC): { double orr, oi; if (((kfn_c)c->imm.p)(RA.r, 0.0, &orr, &oi)) RD.z = orr + oi * I; else RD.z = NAN + NAN * I; } NEXT();
            OP(KERN_CC): { double orr, oi; if (((kfn_c)c->imm.p)(creal(RA.z), cimag(RA.z), &orr, &oi)) RD.z = orr + oi * I; else RD.z = NAN + NAN * I; } NEXT();
            OP(KERN_CR): { double orr, oi; RD.r = ((kfn_c)c->imm.p)(creal(RA.z), cimag(RA.z), &orr, &oi) ? orr : NAN; } NEXT();
            OP(KERN2_RR): { double orr, oi; RD.r = ((kfn_c2)c->imm.p)(RA.r, 0.0, RB.r, 0.0, &orr, &oi) ? orr : NAN; } NEXT();
            OP(KERN2_RC): { double orr, oi; if (((kfn_c2)c->imm.p)(RA.r, 0.0, RB.r, 0.0, &orr, &oi)) RD.z = orr + oi * I; else RD.z = NAN + NAN * I; } NEXT();
            OP(KERN2_CC): { double orr, oi; if (((kfn_c2)c->imm.p)(creal(RA.z), cimag(RA.z), creal(RB.z), cimag(RB.z), &orr, &oi)) RD.z = orr + oi * I; else RD.z = NAN + NAN * I; } NEXT();
            OP(LT_I): RD.i = RA.i < RB.i; NEXT();  OP(LT_R): RD.i = RA.r < RB.r; NEXT();
            OP(LE_I): RD.i = RA.i <= RB.i; NEXT(); OP(LE_R): RD.i = RA.r <= RB.r; NEXT();
            OP(GT_I): RD.i = RA.i > RB.i; NEXT();  OP(GT_R): RD.i = RA.r > RB.r; NEXT();
            OP(GE_I): RD.i = RA.i >= RB.i; NEXT(); OP(GE_R): RD.i = RA.r >= RB.r; NEXT();
            OP(EQ_I): RD.i = RA.i == RB.i; NEXT(); OP(EQ_R): RD.i = RA.r == RB.r; NEXT();
            OP(EQ_C): RD.i = RA.z == RB.z; NEXT();
            OP(NE_I): RD.i = RA.i != RB.i; NEXT(); OP(NE_R): RD.i = RA.r != RB.r; NEXT();
            OP(NE_C): RD.i = RA.z != RB.z; NEXT();
            /* SameQ on machine numbers, matching expr_eq (src/expr.c:622):
             * two NaNs ARE the same, which is what lets an iteration whose
             * orbit reaches NaN terminate instead of spinning.  Complex is
             * componentwise because expr_eq recurses into Complex[re, im]. */
            OP(SAMEQ_R): RD.i = (RA.r == RB.r) || (isnan(RA.r) && isnan(RB.r)); NEXT();
            OP(SAMEQ_C): {
                double ar = creal(RA.z), ai = cimag(RA.z);
                double br = creal(RB.z), bi = cimag(RB.z);
                RD.i = ((ar == br) || (isnan(ar) && isnan(br)))
                    && ((ai == bi) || (isnan(ai) && isnan(bi)));
            } NEXT();
            OP(FAIL): goto vm_fail;
            OP(AND): RD.i = RA.i && RB.i; NEXT();
            OP(OR):  RD.i = RA.i || RB.i; NEXT();
            OP(XOR): RD.i = (!!RA.i) ^ (!!RB.i); NEXT();
            OP(NOT): RD.i = !RA.i; NEXT();
            /* Arbitrary-precision managed scalars.  R[dst].p / R[a].p / R[b].p
             * are container POINTERS bound at frame entry; the op writes THROUGH
             * R[dst].p (never reassigns it) so the binding survives.  The real and
             * complex opcodes need MPFR; without it they can only appear in a
             * program that was never built (the option degrades), so they abort. */
            OP(MG_RCONST):
#ifdef USE_MPFR
                mpfr_set((mpfr_ptr)RD.p, (mpfr_srcptr)c->imm.p, MPFR_RNDN);
#else
                goto vm_fail;
#endif
                NEXT();
            OP(MG_RMOV):
#ifdef USE_MPFR
                mpfr_set((mpfr_ptr)RD.p, (mpfr_srcptr)RA.p, MPFR_RNDN);
#else
                goto vm_fail;
#endif
                NEXT();
            OP(MG_RFROM_I):
#ifdef USE_MPFR
                mpfr_set_si((mpfr_ptr)RD.p, (long)RA.i, MPFR_RNDN);
#else
                goto vm_fail;
#endif
                NEXT();
            OP(MG_RFROM_D):
#ifdef USE_MPFR
                mpfr_set_d((mpfr_ptr)RD.p, RA.r, MPFR_RNDN);
#else
                goto vm_fail;
#endif
                NEXT();
            OP(MG_RUN):
#ifdef USE_MPFR
                ((int(*)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t))c->imm.p)
                    ((mpfr_ptr)RD.p, (mpfr_srcptr)RA.p, MPFR_RNDN);
#else
                goto vm_fail;
#endif
                NEXT();
            OP(MG_RBIN):
#ifdef USE_MPFR
                ((int(*)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t))c->imm.p)
                    ((mpfr_ptr)RD.p, (mpfr_srcptr)RA.p, (mpfr_srcptr)RB.p, MPFR_RNDN);
#else
                goto vm_fail;
#endif
                NEXT();
            OP(MG_RPOWI):
#ifdef USE_MPFR
                mpfr_pow_si((mpfr_ptr)RD.p, (mpfr_srcptr)RA.p, (long)c->imm.i, MPFR_RNDN);
#else
                goto vm_fail;
#endif
                NEXT();
            OP(MG_RCMP):
#ifdef USE_MPFR
                {
                    mpfr_srcptr x = (mpfr_srcptr)RA.p, y = (mpfr_srcptr)RB.p;
                    switch ((long)c->imm.i) {
                        case 0: RD.i = mpfr_less_p(x, y);         break;
                        case 1: RD.i = mpfr_lessequal_p(x, y);    break;
                        case 2: RD.i = mpfr_greater_p(x, y);      break;
                        case 3: RD.i = mpfr_greaterequal_p(x, y); break;
                        case 4: RD.i = mpfr_equal_p(x, y);        break;
                        default: RD.i = !mpfr_equal_p(x, y);      break;
                    }
                }
#else
                goto vm_fail;
#endif
                NEXT();
            OP(MG_ZCONST): mpz_set((mpz_ptr)RD.p, (mpz_srcptr)c->imm.p); NEXT();
            OP(MG_ZMOV):   mpz_set((mpz_ptr)RD.p, (mpz_srcptr)RA.p); NEXT();
            OP(MG_ZFROM_I): mpz_set_si((mpz_ptr)RD.p, (long)RA.i); NEXT();
            OP(MG_ZADD): mpz_add((mpz_ptr)RD.p, (mpz_srcptr)RA.p, (mpz_srcptr)RB.p); NEXT();
            OP(MG_ZSUB): mpz_sub((mpz_ptr)RD.p, (mpz_srcptr)RA.p, (mpz_srcptr)RB.p); NEXT();
            OP(MG_ZMUL): mpz_mul((mpz_ptr)RD.p, (mpz_srcptr)RA.p, (mpz_srcptr)RB.p); NEXT();
            OP(MG_ZNEG): mpz_neg((mpz_ptr)RD.p, (mpz_srcptr)RA.p); NEXT();
            OP(MG_ZPOWI): mpz_pow_ui((mpz_ptr)RD.p, (mpz_srcptr)RA.p, (unsigned long)c->imm.i); NEXT();
            OP(MG_ZCMP): {
                int cmp = mpz_cmp((mpz_srcptr)RA.p, (mpz_srcptr)RB.p);
                switch ((long)c->imm.i) {
                    case 0: RD.i = cmp <  0; break; case 1: RD.i = cmp <= 0; break;
                    case 2: RD.i = cmp >  0; break; case 3: RD.i = cmp >= 0; break;
                    case 4: RD.i = cmp == 0; break; default: RD.i = cmp != 0; break;
                }
            } NEXT();
            OP(MG_CCONST):
#ifdef USE_MPFR
                ncpx_set((ncpx*)RD.p, (const ncpx*)c->imm.p);
#else
                goto vm_fail;
#endif
                NEXT();
            OP(MG_CMOV):
#ifdef USE_MPFR
                ncpx_set((ncpx*)RD.p, (const ncpx*)RA.p);
#else
                goto vm_fail;
#endif
                NEXT();
            OP(MG_CFROM_R):
#ifdef USE_MPFR
                {
                    ncpx* d = (ncpx*)RD.p;
                    mpfr_set(d->re, (mpfr_srcptr)RA.p, MPFR_RNDN);
                    mpfr_set_ui(d->im, 0, MPFR_RNDN);
                }
#else
                goto vm_fail;
#endif
                NEXT();
            OP(MG_CUN):
#ifdef USE_MPFR
                {
                    const ncpx* x = (const ncpx*)RA.p;
                    long sel = (long)c->imm.i;
                    if (sel == MGC_ABS)      ncpx_abs((mpfr_ptr)RD.p, x);
                    else if (sel == MGC_ARG) ncpx_arg((mpfr_ptr)RD.p, x);
                    else {
                        ncpx* d = (ncpx*)RD.p; mpfr_prec_t wp = mpfr_get_prec(d->re);
                        switch (sel) {
                            case MGC_NEG:  ncpx_neg(d, x);       break;
                            case MGC_EXP:  ncpx_exp(d, x, wp);   break;
                            case MGC_LOG:  ncpx_log(d, x, wp);   break;
                            case MGC_SIN:  ncpx_sin(d, x, wp);   break;
                            case MGC_COS:  ncpx_cos(d, x, wp);   break;
                            case MGC_SQRT: ncpx_sqrt(d, x, wp);  break;
                            default: goto vm_fail;
                        }
                    }
                }
#else
                goto vm_fail;
#endif
                NEXT();
            OP(MG_CBIN):
#ifdef USE_MPFR
                {
                    ncpx* d = (ncpx*)RD.p;
                    const ncpx* x = (const ncpx*)RA.p; const ncpx* y = (const ncpx*)RB.p;
                    mpfr_prec_t wp = mpfr_get_prec(d->re);
                    switch ((long)c->imm.i) {
                        case MGC_ADD: ncpx_add(d, x, y);     break;
                        case MGC_SUB: ncpx_sub(d, x, y);     break;
                        case MGC_MUL: ncpx_mul(d, x, y, wp); break;
                        case MGC_DIV: ncpx_div(d, x, y, wp); break;
                        case MGC_POW: ncpx_pow(d, x, y, wp); break;
                        default: goto vm_fail;
                    }
                }
#else
                goto vm_fail;
#endif
                NEXT();
            OP(MG_CPOWI):
#ifdef USE_MPFR
                {
                    ncpx* d = (ncpx*)RD.p; const ncpx* x = (const ncpx*)RA.p;
                    ncpx_pow_d(d, x, (double)(long)c->imm.i, mpfr_get_prec(d->re));
                }
#else
                goto vm_fail;
#endif
                NEXT();
            /* Array ops are out of line: they allocate, they can fail, and
             * keeping them out of the scalar cases costs the scalar path
             * nothing. */
            OP(ARR_FREE): ARROP();
            OP(V_EW):     ARROP();
            OP(V_POW):    ARROP();
            OP(V_KERN):   ARROP();
            OP(V_KERN2):  ARROP();
            OP(V_TOTAL):  ARROP();
            OP(V_NDRED):  ARROP();
            OP(V_NDREDN): ARROP();
            OP(A_NDFN2):  ARROP();
            OP(V_NDFN2):  ARROP();
            OP(V_LEN):    ARROP();
            /* ---- fused elementwise loop (M3b) ----------------------------
             * A_LOAD/A_STORE are the whole point of fusion: an elementwise
             * chain becomes ONE pass over the buffers driven by ordinary scalar
             * bytecode, instead of one full-length ND pass and one temporary
             * buffer per operation.  They are inline (not routed through
             * vm_array_op) because they run once per element.
             *
             * The source dtype is checked per element rather than hoisted: it is
             * a perfectly predicted branch next to the arithmetic it feeds, and
             * hoisting it would mean either specialising the loop at runtime or
             * refusing float32 arrays that work today. */
            OP(A_LOAD_R): {
                const NDArrayData* A_ = &RA.arr->data.ndarray;
                size_t k_ = (size_t)RB.i;
                if (A_->dtype == NDT_FLOAT64)      RD.r = ((const double*)A_->data)[k_];
                else if (A_->dtype == NDT_FLOAT32) RD.r = (double)((const float*)A_->data)[k_];
                else goto vm_fail;   /* promised real, buffer is complex */
            } NEXT();
            OP(A_LOAD_C): {
                const NDArrayData* A_ = &RA.arr->data.ndarray;
                size_t k_ = (size_t)RB.i;
                switch (A_->dtype) {
                    case NDT_FLOAT64:   RD.z = ((const double*)A_->data)[k_]; break;
                    case NDT_FLOAT32:   RD.z = (double)((const float*)A_->data)[k_]; break;
                    case NDT_COMPLEX64: RD.z = ((const double*)A_->data)[2*k_]
                                             + ((const double*)A_->data)[2*k_+1] * I; break;
                    default:            RD.z = (double)((const float*)A_->data)[2*k_]
                                             + (double)((const float*)A_->data)[2*k_+1] * I; break;
                }
            } NEXT();
            /* Integer element access.  Unlike the real form there is no width
             * variation to absorb — a program's integer buffers are always
             * NDT_INT64 (nd_own_copy refuses to build one from anything else),
             * so a foreign dtype here means the promised element type was not
             * kept and the call goes back to the interpreter. */
            OP(A_LOAD_I): {
                const NDArrayData* A_ = &RA.arr->data.ndarray;
                if (A_->dtype != NDT_INT64) goto vm_fail;
                RD.i = ((const int64_t*)A_->data)[(size_t)RB.i];
            } NEXT();
            /* Boolean element access.  One byte per element, held as 0/1 in the
             * integer slot (a compiled boolean is a long long, like True/False).
             * A program's bool buffers are always NDT_BOOL, so a foreign dtype
             * means the promised element type was not kept — back to the
             * interpreter, exactly as A_LOAD_I does for a non-int64 buffer. */
            OP(A_LOAD_B): {
                const NDArrayData* A_ = &RA.arr->data.ndarray;
                if (A_->dtype != NDT_BOOL) goto vm_fail;
                RD.i = ((const uint8_t*)A_->data)[(size_t)RB.i] ? 1 : 0;
            } NEXT();
            OP(A_STORE_I): {
                NDArrayData* A_ = &RD.arr->data.ndarray;
                ((int64_t*)A_->data)[(size_t)RA.i] = RB.i;
            } NEXT();
            OP(A_STORE_B): {
                NDArrayData* A_ = &RD.arr->data.ndarray;
                ((uint8_t*)A_->data)[(size_t)RA.i] = RB.i ? 1 : 0;
            } NEXT();
            OP(A_STORE_R): {
                NDArrayData* A_ = &RD.arr->data.ndarray;
                ((double*)A_->data)[(size_t)RA.i] = RB.r;
            } NEXT();
            OP(A_STORE_C): {
                NDArrayData* A_ = &RD.arr->data.ndarray;
                size_t k_ = (size_t)RA.i;
                ((double*)A_->data)[2*k_]   = creal(RB.z);
                ((double*)A_->data)[2*k_+1] = cimag(RB.z);
            } NEXT();
            /* Append at R[a], growing the buffer when it runs out.  The array
             * was allocated by this call and nothing else holds it, so growing
             * it in place is safe; `dims[0]` is the CAPACITY until A_TRUNC sets
             * the real length. */
            OP(A_PUSH): {
                Expr* A_ = RD.arr;
                if (!A_ || A_->type != EXPR_NDARRAY || A_->data.ndarray.rank != 1) goto vm_fail;
                long long k_ = RA.i;
                if (k_ < 0) goto vm_fail;
                NDArrayData* nd_ = &A_->data.ndarray;
                if (k_ >= nd_->dims[0]) {
                    int64_t cap_ = nd_->dims[0] > 0 ? nd_->dims[0] * 2 : 16;
                    if (cap_ <= k_) cap_ = k_ + 1;
                    if (cap_ > (int64_t)VM_ITER_SAFETY_CAP + 1) goto vm_fail;
                    size_t es_ = ndt_elem_size(nd_->dtype);
                    void* nb_ = realloc(nd_->data, (size_t)cap_ * es_);
                    if (!nb_) goto vm_fail;
                    nd_->data = nb_;
                    nd_->dims[0] = cap_;
                }
                if (AF_R(c->flags) == (unsigned)CT_COMPLEX) {
                    ((double*)nd_->data)[2*(size_t)k_]     = creal(RB.z);
                    ((double*)nd_->data)[2*(size_t)k_ + 1] = cimag(RB.z);
                } else ((double*)nd_->data)[(size_t)k_] = RB.r;
            } NEXT();
            /* Final length, so a capacity the run happened to reach is never
             * visible.  Shrinking in place: a failed realloc keeps the buffer,
             * which is still large enough, so only dims[0] has to be right. */
            OP(A_TRUNC): {
                Expr* A_ = RD.arr;
                if (!A_ || A_->type != EXPR_NDARRAY || A_->data.ndarray.rank != 1) goto vm_fail;
                long long n_ = RA.i;
                NDArrayData* nd_ = &A_->data.ndarray;
                if (n_ < 0 || n_ > nd_->dims[0]) goto vm_fail;
                size_t es_ = ndt_elem_size(nd_->dtype);
                void* nb_ = realloc(nd_->data, (size_t)(n_ ? n_ : 1) * es_);
                if (nb_) nd_->data = nb_;
                nd_->dims[0] = n_;
            } NEXT();
            OP(A_SIZE):     ARROP();
            OP(A_NEWLIKE):  ARROP();
            OP(A_SHAPECHK): ARROP();
            OP(A_COPY):     ARROP();
            OP(A_XFER):     ARROP();
            /* ---- indexed Part (M3c) --------------------------------------
             * One instruction per axis, resolving the subscript against that
             * axis and folding it into the running flat index.  The range check
             * has to be per axis: u[[1, n + 5]] on an n x n array is inside the
             * buffer and reads the row below, which is the one indexing bug a
             * flat bounds check cannot catch. */
            OP(A_AXIS): {
                const Expr* A_ = RB.arr;
                long long ax = c->imm.i;
                if (!A_ || A_->type != EXPR_NDARRAY
                    || ax >= (long long)A_->data.ndarray.rank) goto vm_fail;
                long long len = (long long)A_->data.ndarray.dims[ax];
                long long k_ = RA.i;
                if (k_ < 0) k_ = len + k_ + 1;          /* Part counts from the end */
                if (k_ < 1 || k_ > len) goto vm_fail;   /* -> interpreter, which reports it */
                RD.i = RD.i * len + (k_ - 1);
            } NEXT();
            OP(A_NEW):     do { if (!vm_range_array_op(c, R)) goto vm_fail; } while (0); NEXT();
            OP(A_PART):    do { if (!vm_range_array_op(c, R)) goto vm_fail; } while (0); NEXT();
            OP(A_PARTSET): do { if (!vm_range_array_op(c, R)) goto vm_fail; } while (0); NEXT();
            OP(A_NDFN):    do { if (!vm_range_array_op(c, R)) goto vm_fail; } while (0); NEXT();

            /* ---- Association read ops (B1) --------------------------------
             * `imm.p` is a program-owned AssocSpec; `flags` carries the result
             * element type (LOOKUP) / the KeyFreeQ negation (HASKEY) / the packed
             * dtype (VALUES).  The bag is the spec's constant association, or the
             * borrowed handle in the operand register R[c->a].  A value that does
             * not fit the declared type, or an absent key with no default,
             * declines to the interpreter (goto vm_fail) — the faithful degrade. */
            OP(ASSOC_LOOKUP): {
                const AssocSpec* sp = (const AssocSpec*)c->imm.p;
                Expr* assoc = sp->assoc ? sp->assoc : RA.arr;
                if (!assoc) goto vm_fail;
                Expr* v = assoc_lookup_value(assoc, sp->key);
                if (!v) v = sp->deflt;
                if (!v) goto vm_fail;                      /* absent, no default */
                switch (c->flags) {
                    case CT_INT:  { long long x; if (!cf_to_ll(v, &x)) goto vm_fail; RD.i = x; } break;
                    case CT_REAL: { double x;    if (!cf_to_double(v, &x)) goto vm_fail; RD.r = x; } break;
                    case CT_COMPLEX: { double re, im; if (!cf_to_complex(v, &re, &im)) goto vm_fail;
                                       RD.z = re + im * I; } break;
                    default: goto vm_fail;
                }
            } NEXT();
            OP(ASSOC_HASKEY): {
                const AssocSpec* sp = (const AssocSpec*)c->imm.p;
                Expr* assoc = sp->assoc ? sp->assoc : RA.arr;
                if (!assoc) goto vm_fail;
                long long present = assoc_lookup_value(assoc, sp->key) != NULL;
                RD.i = c->flags ? !present : present;      /* flags=1 => KeyFreeQ */
            } NEXT();
            OP(ASSOC_LEN): {
                Expr* assoc = RA.arr;
                if (!assoc) goto vm_fail;
                RD.i = (long long)assoc->data.function.arg_count;
            } NEXT();
            OP(ASSOC_VALUES): do { if (!vm_assoc_values(c, R)) goto vm_fail; } while (0); NEXT();
            OP(ASSOC_KEYSEL): do { if (!vm_assoc_keysel(c, R)) goto vm_fail; } while (0); NEXT();
            OP(ASSOC_COUNTS): do { if (!vm_assoc_counts(c, R)) goto vm_fail; } while (0); NEXT();
            OP(ASSOC_MAP):    do { if (!vm_assoc_higher(c, R, false)) goto vm_fail; } while (0); NEXT();
            OP(ASSOC_SELECT): do { if (!vm_assoc_higher(c, R, true))  goto vm_fail; } while (0); NEXT();
            OP(ASSOC_SET):    do { if (!vm_assoc_set(c, R)) goto vm_fail; } while (0); NEXT();
            OP(ASSOC_LOOKUP_DYN): {          /* B2: runtime int/real key in R[b] */
                const AssocSpec* sp = (const AssocSpec*)c->imm.p;
                Expr* assoc = sp->assoc ? sp->assoc : RA.arr;
                if (!assoc) goto vm_fail;
                Expr* v = ((unsigned)(c->flags >> 4) & 0xF) == (unsigned)CT_INT
                        ? assoc_lookup_value_i64(assoc, RB.i)
                        : assoc_lookup_value_real(assoc, RB.r);
                if (!v) v = sp->deflt;
                if (!v) goto vm_fail;
                switch ((unsigned)c->flags & 0xF) {
                    case CT_INT:  { long long x; if (!cf_to_ll(v, &x)) goto vm_fail; RD.i = x; } break;
                    case CT_REAL: { double x;    if (!cf_to_double(v, &x)) goto vm_fail; RD.r = x; } break;
                    case CT_COMPLEX: { double re, im; if (!cf_to_complex(v, &re, &im)) goto vm_fail;
                                       RD.z = re + im * I; } break;
                    default: goto vm_fail;
                }
            } NEXT();

            /* ---- strip-mined tile ops (M5b) -------------------------------
             * One opcode, VBLOCK elements, in a loop shaped so the C compiler
             * can vectorise it.  Arithmetic always covers the FULL tile: the
             * load pads a short tail with 1.0, so no op reads uninitialised
             * memory and only load/store/reduce need the live length. */
            #define TD_R  ((double*)RD.p)
            #define TA_R  ((const double*)RA.p)
            #define TB_R  ((const double*)RB.p)
            #define TD_C  ((double _Complex*)RD.p)
            #define TA_C  ((const double _Complex*)RA.p)
            #define TB_C  ((const double _Complex*)RB.p)
            #define VBIN(NAME, TAG, CTY, EXPR) \
                OP(NAME): { CTY* d_ = TD_##TAG; const CTY* a_ = TA_##TAG; \
                            const CTY* b_ = TB_##TAG; (void)a_; (void)b_; \
                            for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = (EXPR); } NEXT()
            #define VUN(NAME, TAG, CTY, EXPR) \
                OP(NAME): { CTY* d_ = TD_##TAG; const CTY* a_ = TA_##TAG; \
                            for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = (EXPR); } NEXT()

            OP(VSETLEN): {                    /* live elements of this tile */
                long long rem = RB.i - RA.i;
                vlen = (int)(rem < VBLOCK ? rem : VBLOCK);
                if (vlen < 0) vlen = 0;
            } NEXT();

            OP(VLOAD_R): {
                const NDArrayData* A_ = &RA.arr->data.ndarray;
                size_t off = (size_t)RB.i;
                double* d_ = TD_R;
                if (A_->dtype == NDT_FLOAT64) {
                    const double* s_ = (const double*)A_->data + off;
                    for (int k_ = 0; k_ < vlen; k_++) d_[k_] = s_[k_];
                } else if (A_->dtype == NDT_FLOAT32) {
                    const float* s_ = (const float*)A_->data + off;
                    for (int k_ = 0; k_ < vlen; k_++) d_[k_] = (double)s_[k_];
                } else goto vm_fail;          /* promised real, buffer is complex */
                for (int k_ = vlen; k_ < VBLOCK; k_++) d_[k_] = 1.0;   /* safe pad */
            } NEXT();
            OP(VLOAD_C): {
                const NDArrayData* A_ = &RA.arr->data.ndarray;
                size_t off = (size_t)RB.i;
                double _Complex* d_ = TD_C;
                switch (A_->dtype) {
                    case NDT_FLOAT64: { const double* s_ = (const double*)A_->data + off;
                        for (int k_ = 0; k_ < vlen; k_++) d_[k_] = s_[k_];
                        break; }
                    case NDT_FLOAT32: { const float* s_ = (const float*)A_->data + off;
                        for (int k_ = 0; k_ < vlen; k_++) d_[k_] = (double)s_[k_];
                        break; }
                    case NDT_COMPLEX64: { const double* s_ = (const double*)A_->data + 2 * off;
                        for (int k_ = 0; k_ < vlen; k_++) d_[k_] = s_[2*k_] + s_[2*k_+1] * I;
                        break; }
                    default: { const float* s_ = (const float*)A_->data + 2 * off;
                        for (int k_ = 0; k_ < vlen; k_++)
                            d_[k_] = (double)s_[2*k_] + (double)s_[2*k_+1] * I;
                        break; }
                }
                for (int k_ = vlen; k_ < VBLOCK; k_++) d_[k_] = 1.0;
            } NEXT();
            /* The store is where a fused loop honours the element-type promise.
             * A real-typed program that produced a non-finite element is exactly
             * the case where the interpreter would have gone complex (Sqrt of a
             * negative, Log of a negative) or hit a pole, and the delegated path
             * fails there too because its kernels reject non-finite elements.
             * Returning a buffer full of NaN instead would silently hand back a
             * real array where the interpreter gives a complex one.
             *
             * The check is accumulated branchlessly and tested once per tile, so
             * it does not break the vectorisation of the copy. */
            OP(VSTORE_R): {
                double* d_ = (double*)RD.arr->data.ndarray.data + (size_t)RA.i;
                const double* s_ = TB_R;
                int bad_ = 0;
                for (int k_ = 0; k_ < vlen; k_++) {
                    double x_ = s_[k_];
                    d_[k_] = x_;
                    bad_ |= !(x_ - x_ == 0.0);          /* NaN or +-Inf */
                }
                if (bad_) goto vm_fail;
            } NEXT();
            OP(VSTORE_C): {
                double* d_ = (double*)RD.arr->data.ndarray.data + 2 * (size_t)RA.i;
                const double _Complex* s_ = TB_C;
                int bad_ = 0;
                for (int k_ = 0; k_ < vlen; k_++) {
                    double re_ = creal(s_[k_]), im_ = cimag(s_[k_]);
                    d_[2*k_] = re_; d_[2*k_+1] = im_;
                    bad_ |= !(re_ - re_ == 0.0) | !(im_ - im_ == 0.0);
                }
                if (bad_) goto vm_fail;
            } NEXT();
            /* Sequential accumulation, deliberately: it reproduces the
             * element-at-a-time fused loop's summation order exactly, so
             * strip-mining does not change any already-tested result. */
            OP(VSUM_R): { const double* s_ = TA_R; double acc_ = RD.r;
                          for (int k_ = 0; k_ < vlen; k_++) acc_ += s_[k_];
                          RD.r = acc_; } NEXT();
            OP(VSUM_C): { const double _Complex* s_ = TA_C; double _Complex acc_ = RD.z;
                          for (int k_ = 0; k_ < vlen; k_++) acc_ += s_[k_];
                          RD.z = acc_; } NEXT();

            OP(VSPLAT_R): { double* d_ = TD_R; double v_ = RA.r;
                            for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = v_; } NEXT();
            OP(VSPLAT_C): { double _Complex* d_ = TD_C; double _Complex v_ = RA.z;
                            for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = v_; } NEXT();
            OP(VR2C): { double _Complex* d_ = TD_C; const double* a_ = TA_R;
                        for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = a_[k_]; } NEXT();

            VBIN(VADD_R, R, double, a_[k_] + b_[k_]);
            VBIN(VSUB_R, R, double, a_[k_] - b_[k_]);
            VBIN(VMUL_R, R, double, a_[k_] * b_[k_]);
            VBIN(VDIV_R, R, double, a_[k_] / b_[k_]);
            VUN(VNEG_R, R, double, -a_[k_]);
            VUN(VINV_R, R, double, 1.0 / a_[k_]);
            VBIN(VADD_C, C, double _Complex, a_[k_] + b_[k_]);
            VBIN(VSUB_C, C, double _Complex, a_[k_] - b_[k_]);
            VBIN(VMUL_C, C, double _Complex, a_[k_] * b_[k_]);
            VBIN(VDIV_C, C, double _Complex, a_[k_] / b_[k_]);
            VUN(VNEG_C, C, double _Complex, -a_[k_]);
            VUN(VINV_C, C, double _Complex, 1.0 / a_[k_]);
            /* Small exponents are unrolled with EXACTLY the association
             * ipow_* uses (x^2 = x*x, x^3 = x*(x*x), x^4 = (x*x)*(x*x)), so the
             * result is bitwise identical to the general path while the tile
             * loop becomes a vectorisable multiply instead of a per-element
             * function call with a loop inside it. */
            OP(VPOWI_R): { double* d_ = TD_R; const double* a_ = TA_R; long long e_ = c->imm.i;
                if (e_ == 2)      for (int k_ = 0; k_ < VBLOCK; k_++) { double t_ = a_[k_]; d_[k_] = t_ * t_; }
                else if (e_ == 3) for (int k_ = 0; k_ < VBLOCK; k_++) { double t_ = a_[k_]; d_[k_] = t_ * (t_ * t_); }
                else if (e_ == 4) for (int k_ = 0; k_ < VBLOCK; k_++) { double t_ = a_[k_], s_ = t_ * t_; d_[k_] = s_ * s_; }
                else              for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = ipow_r(a_[k_], e_);
            } NEXT();
            OP(VPOWI_C): { double _Complex* d_ = TD_C; const double _Complex* a_ = TA_C; long long e_ = c->imm.i;
                if (e_ == 2)      for (int k_ = 0; k_ < VBLOCK; k_++) { double _Complex t_ = a_[k_]; d_[k_] = t_ * t_; }
                else if (e_ == 3) for (int k_ = 0; k_ < VBLOCK; k_++) { double _Complex t_ = a_[k_]; d_[k_] = t_ * (t_ * t_); }
                else if (e_ == 4) for (int k_ = 0; k_ < VBLOCK; k_++) { double _Complex t_ = a_[k_], s_ = t_ * t_; d_[k_] = s_ * s_; }
                else              for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = ipow_c(a_[k_], e_);
            } NEXT();
            VBIN(VPOW_R, R, double, pow(a_[k_], b_[k_]));
            VBIN(VPOW_C, C, double _Complex, cpow(a_[k_], b_[k_]));
            VBIN(VATAN2_R, R, double, atan2(b_[k_], a_[k_]));
            VUN(VSQRT_R, R, double, sqrt(a_[k_]));   VUN(VSQRT_C, C, double _Complex, csqrt(a_[k_]));
            VUN(VEXP_R, R, double, exp(a_[k_]));    VUN(VEXP_C, C, double _Complex, cexp(a_[k_]));
            VUN(VLOG_R, R, double, log(a_[k_]));    VUN(VLOG_C, C, double _Complex, clog(a_[k_]));
            VUN(VSIN_R, R, double, sin(a_[k_]));    VUN(VSIN_C, C, double _Complex, csin(a_[k_]));
            VUN(VCOS_R, R, double, cos(a_[k_]));    VUN(VCOS_C, C, double _Complex, ccos(a_[k_]));
            VUN(VTAN_R, R, double, tan(a_[k_]));    VUN(VTAN_C, C, double _Complex, ctan(a_[k_]));
            VUN(VSINH_R, R, double, sinh(a_[k_]));   VUN(VSINH_C, C, double _Complex, csinh(a_[k_]));
            VUN(VCOSH_R, R, double, cosh(a_[k_]));   VUN(VCOSH_C, C, double _Complex, ccosh(a_[k_]));
            VUN(VTANH_R, R, double, tanh(a_[k_]));   VUN(VTANH_C, C, double _Complex, ctanh(a_[k_]));
            VUN(VASIN_R, R, double, asin(a_[k_]));   VUN(VASIN_C, C, double _Complex, casin(a_[k_]));
            VUN(VACOS_R, R, double, acos(a_[k_]));   VUN(VACOS_C, C, double _Complex, cacos(a_[k_]));
            VUN(VATAN_R, R, double, atan(a_[k_]));   VUN(VATAN_C, C, double _Complex, catan(a_[k_]));
            VUN(VABS_R, R, double, fabs(a_[k_]));
            VUN(VSIGN_R, R, double, (double)((a_[k_] > 0) - (a_[k_] < 0)));
            VUN(VERF_R, R, double, erf(a_[k_]));
            VUN(VERFC_R, R, double, erfc(a_[k_]));
            /* complex tile in, real tile out */
            OP(VABS_C): { double* d_ = TD_R; const double _Complex* a_ = TA_C;
                          for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = cabs(a_[k_]); } NEXT();
            OP(VRE_C):  { double* d_ = TD_R; const double _Complex* a_ = TA_C;
                          for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = creal(a_[k_]); } NEXT();
            OP(VIM_C):  { double* d_ = TD_R; const double _Complex* a_ = TA_C;
                          for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = cimag(a_[k_]); } NEXT();
            OP(VARG_C): { double* d_ = TD_R; const double _Complex* a_ = TA_C;
                          for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = carg(a_[k_]); } NEXT();
            VUN(VCONJ_C, C, double _Complex, conj(a_[k_]));

            /* Kernel tiles: the indirect call per element blocks vectorisation,
             * but these are libm-class kernels where the call dominates anyway —
             * the win here is purely the amortised dispatch. */
            OP(VKERN_RR): { double* d_ = TD_R; const double* a_ = TA_R; double o_;
                for (int k_ = 0; k_ < VBLOCK; k_++)
                    d_[k_] = ((kfn_r)c->imm.p)(a_[k_], &o_) ? o_ : NAN; } NEXT();
            OP(VKERN_R2R): { double* d_ = TD_R; const double* a_ = TA_R; double or_, oi_;
                for (int k_ = 0; k_ < VBLOCK; k_++)
                    d_[k_] = ((kfn_c)c->imm.p)(a_[k_], 0.0, &or_, &oi_) ? or_ : NAN; } NEXT();
            OP(VKERN_RC): { double _Complex* d_ = TD_C; const double* a_ = TA_R; double or_, oi_;
                for (int k_ = 0; k_ < VBLOCK; k_++)
                    d_[k_] = ((kfn_c)c->imm.p)(a_[k_], 0.0, &or_, &oi_) ? or_ + oi_ * I : NAN + NAN * I; } NEXT();
            OP(VKERN_CC): { double _Complex* d_ = TD_C; const double _Complex* a_ = TA_C; double or_, oi_;
                for (int k_ = 0; k_ < VBLOCK; k_++)
                    d_[k_] = ((kfn_c)c->imm.p)(creal(a_[k_]), cimag(a_[k_]), &or_, &oi_)
                             ? or_ + oi_ * I : NAN + NAN * I; } NEXT();
            OP(VKERN_CR): { double* d_ = TD_R; const double _Complex* a_ = TA_C; double or_, oi_;
                for (int k_ = 0; k_ < VBLOCK; k_++)
                    d_[k_] = ((kfn_c)c->imm.p)(creal(a_[k_]), cimag(a_[k_]), &or_, &oi_) ? or_ : NAN; } NEXT();
            OP(VKERN2_RR): { double* d_ = TD_R; const double* a_ = TA_R; const double* b_ = TB_R; double or_, oi_;
                for (int k_ = 0; k_ < VBLOCK; k_++)
                    d_[k_] = ((kfn_c2)c->imm.p)(a_[k_], 0.0, b_[k_], 0.0, &or_, &oi_) ? or_ : NAN; } NEXT();
            OP(VKERN2_RC): { double _Complex* d_ = TD_C; const double* a_ = TA_R; const double* b_ = TB_R; double or_, oi_;
                for (int k_ = 0; k_ < VBLOCK; k_++)
                    d_[k_] = ((kfn_c2)c->imm.p)(a_[k_], 0.0, b_[k_], 0.0, &or_, &oi_)
                             ? or_ + oi_ * I : NAN + NAN * I; } NEXT();
            OP(VKERN2_CC): { double _Complex* d_ = TD_C; const double _Complex* a_ = TA_C;
                             const double _Complex* b_ = TB_C; double or_, oi_;
                for (int k_ = 0; k_ < VBLOCK; k_++)
                    d_[k_] = ((kfn_c2)c->imm.p)(creal(a_[k_]), cimag(a_[k_]), creal(b_[k_]), cimag(b_[k_]),
                                                &or_, &oi_) ? or_ + oi_ * I : NAN + NAN * I; } NEXT();
            #undef VBIN
            #undef VUN
            /* A compiled callee, invoked directly: machine values in, machine
             * value out, no Expr and no evaluator round-trip.  This is what
             * lifts the inliner's depth cap — beyond it the callee is CALLED
             * rather than pasted in, so deep chains compile instead of bailing. */
            OP(CALL): { if (!vm_call((const CompiledProgram*)c->imm.p,
                                     &RA, c->flags, &RD)) goto vm_fail; } NEXT();
            /* An n-ary machine kernel: arguments in `flags` consecutive
             * registers starting at `a`, the kernel pointer in the immediate. */
            OP(KERNN): {
                double ar[8], ai[8], orr, oi;
                unsigned na_ = c->flags;
                if (na_ > 8) goto vm_fail;
                for (unsigned k_ = 0; k_ < na_; k_++) { ar[k_] = R[c->a + k_].r; ai[k_] = 0.0; }
                if (!((kfn_cn)c->imm.p)(ar, ai, (size_t)na_, &orr, &oi)) RD.r = NAN;
                else RD.r = orr;
            } NEXT();
            OP(NOP): NEXT();
            OP(RET): return;
#if !VM_THREADED
            default: return;
        }
        pc++;
    }
    return;
#endif
vm_fail:
    *failed = true;   /* abort: the caller releases live array temporaries */
    return;
    #undef OP
    #undef NEXT
    #undef JUMP
    #undef ARROP
}
long compiled_prec_bits(const CompiledProgram* p) { return p ? p->prec_bits : 0; }

CompileType compiled_result_type(const CompiledProgram* p) { return p->result_type; }
bool        compiled_result_built(const CompiledProgram* p) { return p && p->result_built; }
size_t compiled_num_args(const CompiledProgram* p) { return p->nargs; }
size_t compiled_num_instructions(const CompiledProgram* p) { return p->n; }
bool compiled_program_all_real(const CompiledProgram* p) { return p && p->all_real; }
size_t compiled_num_cse(const CompiledProgram* p) { return p ? (size_t)p->ncse : 0; }

size_t compiled_arg_deps(const CompiledProgram* p, int* deps, size_t cap) {
    size_t n = 0;
    for (size_t k = 0; k < p->nargs && n < cap; k++) if (p->argdep[k]) deps[n++] = (int)k;
    return n;
}

#ifdef USE_MPFR
static bool mgd_ncpx_set_from_expr(ncpx* z, const Expr* e);   /* defined below */
#endif

/* Bind one boxed argument to its register.  Returns false if the value does not
 * match the declared type, so the caller can fall back to the interpreter.
 * An array argument is BORROWED: its handle is stored, never freed. */
static bool load_arg(Slot* s, const CompileValue* v, CompileType want) {
    if (CT_IS_ARRAY(want)) {
        Expr* x = CT_IS_ARRAY(v->type) ? v->v.a : NULL;
        if (!x || x->type != EXPR_NDARRAY || x->data.ndarray.rank != CT_RANK(want))
            return false;
        s->arr = x;
        return true;
    }
    if (CT_IS_ASSOC(want)) {
        /* A read-only Association bag rides in `.arr` of a SCALAR-bank argument
         * register (arg registers are < arr_base), so arr_sweep — which only
         * walks [arr_base, tile_base) — never frees this borrowed handle. */
        Expr* x = CT_IS_ASSOC(v->type) ? v->v.a : NULL;
        if (!x || !is_association(x)) return false;
        s->arr = x;
        return true;
    }
    if (CT_IS_MANAGED(want)) {
        /* The container is already bound in s->p by managed_frame_enter (which
         * runs before load_arg); set its VALUE from the boxed numeric Expr in
         * v.a, rounded to the container's fixed working precision. */
        const Expr* e = v->v.a;
        if (!e) return false;
        if (want == CT_BIGINT) return mgd_mpz_from_expr((mpz_ptr)s->p, e);
#ifdef USE_MPFR
        if (want == CT_BIGREAL) return mgd_mpfr_from_expr((mpfr_ptr)s->p, e);
        if (want == CT_BIGCOMPLEX) return mgd_ncpx_set_from_expr((ncpx*)s->p, e);
#endif
        return false;
    }
    if (CT_IS_ARRAY(v->type) || CT_IS_ASSOC(v->type)) return false;
    /* coerce the boxed arg to the register's declared type */
    switch (v->type) {
        case CT_BOOL: s->i = v->v.b ? 1 : 0; return true;
        case CT_INT:  if (want == CT_INT) { s->i = v->v.i; } else if (want == CT_REAL) { s->r = (double)v->v.i; } else { s->z = (double)v->v.i; } return true;
        case CT_REAL: if (want == CT_COMPLEX) { s->z = v->v.r; } else { s->r = v->v.r; } return true;
        case CT_COMPLEX: s->z = v->v.z; return true;
        default: return false;
    }
}

static bool finite_result(const Slot* s, CompileType t) {
    if (t == CT_REAL) return isfinite(s->r);
    if (t == CT_COMPLEX) return isfinite(creal(s->z)) && isfinite(cimag(s->z));
    return true;   /* int/bool/array always ok */
}

/* Teardown for the array bank.  Every array register is NULLed before a call
 * (the frame is reused, so a handle left over from a previous call must never
 * be touched) and any that is still live afterwards is released here — the
 * belt to OP_ARR_FREE's braces, and the only cleanup on the abort path.
 * A result array is NULLed out by the caller first, so it survives. */
/* ------------------------------------------------------------------ *
 *  Call frames                                                        *
 * ------------------------------------------------------------------ *
 * The frame used to be a single buffer owned by the CompiledProgram and mutated
 * through a `const CompiledProgram*`, which meant a program was neither
 * reentrant nor thread-safe: two live calls of the same program shared one
 * register file, and (once fusion landed) one set of tile buffers.
 *
 * Frames now come from the C stack, which is per-thread and naturally nested, so
 * reentrancy and thread-safety both fall out with no arena to own, grow or free
 * — and no allocation at all for any program that fits the fixed buffer, which
 * is every realistic one.  Larger programs fall back to a single malloc.
 *
 * Layout: `nreg` register Slots, then `ntiles * VBLOCK` Slots of tile storage.
 * A Slot is exactly `sizeof(double _Complex)`, so tile storage is correctly
 * aligned for both real and complex tiles.  VM_STACK_SLOTS is defined up with
 * the parallel strip loop, which needs the same frame shape for its workers. */

/* Point each tile register at its slice of the frame's tile storage. */
static void frame_bind_tiles(const CompiledProgram* p, Slot* R) {
    Slot* tiles = R + p->nreg;
    for (int k = 0; k < p->ntiles; k++)
        R[p->tile_base + k].p = tiles + (size_t)k * VBLOCK;
}

/* Only the ARRAY bank: tile slots hold pointers into the frame's own storage,
 * so clearing them would strand the tiles for the rest of the call. */
static void arr_reset(const CompiledProgram* p, Slot* R) {
    for (int r = p->arr_base; r < p->tile_base; r++) R[r].arr = NULL;
}
static void arr_sweep(const CompiledProgram* p, Slot* R) {
    for (int r = p->arr_base; r < p->tile_base; r++)
        if (R[r].arr) { expr_free(R[r].arr); R[r].arr = NULL; }
}

/* ---- Arbitrary-precision container lifecycle (per call) ----
 * A managed register's Slot holds a pointer to a heap container (mpz_t / mpfr_t
 * / ncpx) that needs init/clear.  These bind them at frame entry and release
 * them at exit — the container analogue of arr_reset/arr_sweep.  Every path is
 * gated on nmgd_slots, so a machine program (nmgd_slots == 0) never enters here:
 * the machine path pays nothing.
 *
 * The containers live in a WARM thread-local arena (design doc §7): frame entry
 * ACQUIRES a run of cells (initialising a cell lazily, or re-initialising it on a
 * type/precision change) and frame exit RETREATS the arena stack pointer WITHOUT
 * clearing, so the containers — and their heap limb buffers — stay allocated for
 * the next call.  After warm-up a per-call cost is O(nmgd) pointer binds, no
 * malloc/init/clear/free, which is what makes the per-sample sampler path fast at
 * high precision.  Per-thread, so reentrancy nests through the stack pointer and
 * worker threads never share; gated on nmgd_slots, so the machine path (which has
 * none) never touches any of it. */
#ifdef USE_MPFR
static bool mgd_ncpx_set_from_expr(ncpx* z, const Expr* e) {
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im))
        return mgd_mpfr_from_expr(z->re, re) && mgd_mpfr_from_expr(z->im, im);
    if (!mgd_mpfr_from_expr(z->re, e)) return false;  /* a real value */
    mpfr_set_ui(z->im, 0, MPFR_RNDN);
    return true;
}
#endif

/* One warm arena cell: a container of a known kind and (for MPFR/ncpx) precision,
 * kept alive across calls. `kind`: -1 unused, 0 mpfr, 1 mpz, 2 ncpx. */
typedef struct { int kind; long prec; void* ptr; } MgdCell;
typedef struct { MgdCell* cells; int cap; int sp; } MgdArena;
static VM_TLS MgdArena g_mgd_arena;

static int mgd_kind_of(CompileType t) {
    return t == CT_BIGINT ? 1 : t == CT_BIGCOMPLEX ? 2 : 0;
}
static void* mgd_cell_alloc(int kind, long prec) {
    (void)prec;
    switch (kind) {
        case 1: { mpz_ptr z = malloc(sizeof *z); if (z) mpz_init(z); return z; }
#ifdef USE_MPFR
        case 0: { mpfr_ptr m = malloc(sizeof *m); if (m) mpfr_init2(m, (mpfr_prec_t)prec); return m; }
        case 2: { ncpx* z = malloc(sizeof *z); if (z) ncpx_init(z, (mpfr_prec_t)prec); return z; }
#endif
        default: return NULL;
    }
}
static void mgd_cell_free(int kind, void* ptr) {
    if (!ptr) return;
    switch (kind) {
        case 1: mpz_clear((mpz_ptr)ptr); break;
#ifdef USE_MPFR
        case 0: mpfr_clear((mpfr_ptr)ptr); break;
        case 2: ncpx_clear((ncpx*)ptr); break;
#endif
        default: break;
    }
    free(ptr);
}
/* Set the precision of an already-allocated MPFR/ncpx cell (kept struct, resized
 * limb buffer) — cheap when the precision is stable, which is the common case. */
static void mgd_cell_reprec(int kind, void* ptr, long prec) {
    (void)kind; (void)ptr; (void)prec;
#ifdef USE_MPFR
    if (kind == 0) mpfr_set_prec((mpfr_ptr)ptr, (mpfr_prec_t)prec);
    else if (kind == 2) {
        ncpx* z = (ncpx*)ptr;
        mpfr_set_prec(z->re, (mpfr_prec_t)prec);
        mpfr_set_prec(z->im, (mpfr_prec_t)prec);
    }
#endif
}

static void managed_frame_exit(const CompiledProgram* p, Slot* R) {
    /* Retreat the arena stack pointer; the containers stay warm for the next
     * call.  Null the frame slots so a stale pointer into the arena is never read
     * (the frame is per-call C-stack storage, so this is belt-and-braces). */
    g_mgd_arena.sp -= p->nmgd_slots;
    if (g_mgd_arena.sp < 0) g_mgd_arena.sp = 0;
    for (int i = 0; i < p->nmgd_slots; i++) R[p->mgd_slots[i].reg].p = NULL;
}
static bool managed_frame_enter(const CompiledProgram* p, Slot* R) {
    MgdArena* a = &g_mgd_arena;
    int base = a->sp;
    int need = base + p->nmgd_slots;
    if (need > a->cap) {
        int nc = a->cap ? a->cap * 2 : 16;
        while (nc < need) nc *= 2;
        MgdCell* nb = realloc(a->cells, (size_t)nc * sizeof *nb);
        if (!nb) return false;
        for (int i = a->cap; i < nc; i++) { nb[i].kind = -1; nb[i].prec = 0; nb[i].ptr = NULL; }
        a->cells = nb; a->cap = nc;
    }
    for (int i = 0; i < p->nmgd_slots; i++) {
        MgdCell* cell = &a->cells[base + i];
        int  want_kind = mgd_kind_of(p->mgd_slots[i].type);
        long want_prec = p->prec_bits;
        if (cell->ptr == NULL) {                     /* first use of this cell */
            cell->ptr = mgd_cell_alloc(want_kind, want_prec);
            if (!cell->ptr) { a->sp = base; return false; }
            cell->kind = want_kind; cell->prec = want_prec;
        } else if (cell->kind != want_kind) {        /* reused as a different type */
            mgd_cell_free(cell->kind, cell->ptr);
            cell->ptr = mgd_cell_alloc(want_kind, want_prec);
            if (!cell->ptr) { cell->kind = -1; a->sp = base; return false; }
            cell->kind = want_kind; cell->prec = want_prec;
        } else if (want_kind != 1 && cell->prec != want_prec) {  /* precision change */
            mgd_cell_reprec(want_kind, cell->ptr, want_prec);
            cell->prec = want_prec;
        }
        R[p->mgd_slots[i].reg].p = cell->ptr;
    }
    a->sp = base + p->nmgd_slots;
    return true;
}

/* Build the result Expr from a managed result container, or NULL on a non-finite
 * value (mirrors the finite_result gate — the caller then interprets). */
static Expr* managed_result_expr(const CompiledProgram* p, const Slot* r) {
    switch ((int)p->result_type) {
        case CT_BIGINT: {
            Expr* e = expr_new_bigint_from_mpz((mpz_srcptr)r->p);
            return expr_bigint_normalize(e);
        }
#ifdef USE_MPFR
        case CT_BIGREAL: {
            mpfr_srcptr m = (mpfr_srcptr)r->p;
            if (!mpfr_number_p(m)) return NULL;   /* nan/inf -> fall back */
            return expr_new_mpfr_copy(m);
        }
        case CT_BIGCOMPLEX: {
            ncpx* z = (ncpx*)r->p;
            if (!mpfr_number_p(z->re) || !mpfr_number_p(z->im)) return NULL;
            return numeric_mpfr_make_complex(z->re, z->im);
        }
#endif
        default: return NULL;
    }
}

/* Depth guard for OP_CALL.  Frames live on the C stack, so unbounded nesting
 * would overflow it rather than fail cleanly; a compiled program that recurses
 * past this simply fails and the caller falls back to the interpreter. */
#define VM_MAX_CALL_DEPTH 200
static VM_TLS int vm_call_depth = 0;

/* Run `cp` on `nargs` machine values taken straight from the caller's registers.
 * The caller coerced them to the callee's declared types at emit time, so the
 * copy is a raw Slot move — no boxing, no Expr, no evaluator.  The callee gets
 * its OWN frame, which is what makes a compiled program reentrant. */
static bool vm_call(const CompiledProgram* cp, const Slot* argv, unsigned nargs, Slot* dst) {
    if (!cp || cp->nargs != (size_t)nargs) return false;
    if (vm_call_depth >= VM_MAX_CALL_DEPTH) return false;
    /* A managed (arbitrary-precision) callee is not reachable in v1 — emit_mgd
     * never emits OP_CALL, and the raw-Slot calling convention here cannot carry
     * container-typed arguments — so decline rather than misread a container. */
    if (cp->nmgd_slots) return false;

    Slot stackframe[VM_STACK_SLOTS];
    Slot* heap = NULL;
    Slot* R = stackframe;
    if (cp->frame_slots > VM_STACK_SLOTS) {
        heap = malloc(cp->frame_slots * sizeof(Slot));
        if (!heap) return false;
        R = heap;
    }
    if (cp->ntiles) frame_bind_tiles(cp, R);
    for (unsigned k = 0; k < nargs; k++) R[k] = argv[k];
    arr_reset(cp, R);

    vm_call_depth++;
    bool failed = false;
    vm_run(cp->code, cp->n, R, &failed);
    vm_call_depth--;

    Slot* r = &R[cp->result_reg];
    bool good = !failed && finite_result(r, cp->result_type)
                && !CT_IS_ARRAY(cp->result_type);
    if (good) *dst = *r;
    arr_sweep(cp, R);
    free(heap);
    return good;
}

bool compiled_eval(const CompiledProgram* p, const CompileValue* args, CompileValue* out) {
    Slot stackframe[VM_STACK_SLOTS];
    Slot* heap = NULL;
    Slot* R = stackframe;
    if (p->frame_slots > VM_STACK_SLOTS) {
        heap = malloc(p->frame_slots * sizeof(Slot));
        if (!heap) return false;
        R = heap;
    }
    if (p->ntiles) frame_bind_tiles(p, R);
    if (p->nmgd_slots && !managed_frame_enter(p, R)) { free(heap); return false; }

    for (size_t k = 0; k < p->nargs; k++)
        if (!load_arg(&R[k], &args[k], p->arg_types[k])) {
            if (p->nmgd_slots) managed_frame_exit(p, R);
            free(heap); return false;
        }
    arr_reset(p, R);
    bool failed = false;
    vm_run(p->code, p->n, R, &failed);
    Slot* r = &R[p->result_reg];
    out->type = p->result_type;
    /* An array OR an owned association result (B3) lives in the .arr slot; both
     * transfer ownership to the caller and both leave the array bank to sweep. */
    if (CT_IS_ARRAY(p->result_type) || CT_IS_ASSOC(p->result_type)) {
        out->v.a = failed ? NULL : r->arr;
        if (!failed) r->arr = NULL;    /* ownership transfers to the caller */
        arr_sweep(p, R);
        if (p->nmgd_slots) managed_frame_exit(p, R);
        free(heap);
        return !failed && out->v.a != NULL;
    }
    /* An arbitrary-precision result is COPIED out of its container into a fresh
     * Expr before the containers are cleared; the caller owns it, as for arrays. */
    if (CT_IS_MANAGED(p->result_type)) {
        Expr* rv = failed ? NULL : managed_result_expr(p, r);
        out->v.a = rv;
        arr_sweep(p, R);
        managed_frame_exit(p, R);
        free(heap);
        return rv != NULL;
    }
    switch (p->result_type) {
        case CT_BOOL: out->v.b = (unsigned char)(r->i != 0); break;
        case CT_INT:  out->v.i = r->i; break;
        case CT_REAL: out->v.r = r->r; break;
        case CT_COMPLEX: out->v.z = r->z; break;
        default: break;
    }
    arr_sweep(p, R);
    bool good = !failed && finite_result(r, p->result_type);
    if (p->nmgd_slots) managed_frame_exit(p, R);
    free(heap);
    return good;
}

bool compiled_eval_real(const CompiledProgram* p, const double* args, double* out) {
    if (!p->all_real) return false;   /* implies no array registers and no tiles */
    Slot stackframe[VM_STACK_SLOTS];
    Slot* heap = NULL;
    Slot* R = stackframe;
    if (p->frame_slots > VM_STACK_SLOTS) {
        heap = malloc(p->frame_slots * sizeof(Slot));
        if (!heap) return false;
        R = heap;
    }
    for (size_t k = 0; k < p->nargs; k++) R[k].r = args[k];
    bool failed = false;
    vm_run(p->code, p->n, R, &failed);
    *out = R[p->result_reg].r;
    free(heap);
    /* `failed` was computed and dropped here for as long as no opcode an
     * all-Real program could contain was able to abort — array ops can, and an
     * all-Real program has no array registers.  OP_FAIL changed that: an
     * unbounded FixedPoint hitting the safety cap must drop the call to the
     * interpreter, not hand back whatever the accumulator happened to hold. */
    return !failed && isfinite(*out);
}

bool compiled_eval_real_batch(const CompiledProgram* const* progs, size_t nprogs,
                              const double* args, size_t nargs, double* out) {
    if (nprogs == 0) return true;
    /* one shared frame = the widest program's (big enough for every program's
     * registers); the argument region [0,nargs) is loaded once and never written
     * by any program (dst registers are always temporaries >= nargs). */
    size_t widest = 0;
    for (size_t i = 0; i < nprogs; i++)
        if (progs[i]->frame_slots > widest) widest = progs[i]->frame_slots;

    Slot stackframe[VM_STACK_SLOTS];
    Slot* heap = NULL;
    Slot* F = stackframe;
    if (widest > VM_STACK_SLOTS) {
        heap = malloc(widest * sizeof(Slot));
        if (!heap) return false;
        F = heap;
    }
    for (size_t k = 0; k < nargs; k++) F[k].r = args[k];
    for (size_t i = 0; i < nprogs; i++) {
        /* all_real implies no array registers and no tiles, so one shared frame
         * is enough and needs no per-program tile binding. */
        if (!progs[i]->all_real) { free(heap); return false; }
        bool failed = false;
        vm_run(progs[i]->code, progs[i]->n, F, &failed);
        out[i] = F[(size_t)progs[i]->result_reg].r;
        if (failed || !isfinite(out[i])) { free(heap); return false; }
    }
    free(heap);
    return true;
}

void compiled_free(CompiledProgram* p) {
    if (!p) return;
    for (int i = 0; i < p->nploops; i++) free(p->ploops[i].code);
    free(p->ploops);
    for (int i = 0; i < p->nparts; i++) compile_partspec_free(p->parts[i]);
    free(p->parts);
    for (int i = 0; i < p->nassocs; i++) compile_assocspec_free(p->assocs[i]);
    free(p->assocs);
    for (int i = 0; i < p->ncallees; i++) compile_assoccallee_free(p->callees[i]);
    free(p->callees);
    for (int i = 0; i < p->nmgd_consts; i++) mgd_const_free(&p->mgd_consts[i]);
    free(p->mgd_consts);
    free(p->mgd_slots);
    free(p->code); free(p->arg_types); free(p->argdep);
    free(p);
}
