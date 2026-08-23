/*
 * reduce_form.h
 *
 * Internal logical normal-form layer for `Reduce` (REDUCE_PLAN.md §"Internal
 * solution normal form").  Reduce's engines all *produce* a disjunction of
 * cell-conjunctions, and a real Reduce result is overwhelmingly a top-level
 * `Or` of `And`s -- so the working representation is DNF:
 *
 *     RForm  = Or  of RConj      (a disjunction; `is_true` short-circuits)
 *     RConj  = And of RAtom      (a conjunction;  `is_false` short-circuits)
 *     RAtom  = (poly) REL 0      (one canonical relational atom)
 *
 * Constant sentinels:
 *     True   -> RForm with is_true == true
 *     False  -> RForm with n == 0 and is_true == false  (an empty disjunction)
 *
 * Greater / GreaterEqual are canonicalised away (operands swapped) so only
 * LT / LE / EQ / NE / ELEM ever reach an engine, halving the case matrix.
 *
 * Ownership: an RAtom owns its `poly` (and, for R_ELEM, `elem_dom`).  Pushing an
 * RAtom into an RConj transfers that ownership.  rform_and / rform_or CONSUME
 * both operands and return a freshly-owned form.  Everything is freed through
 * rform_free -> rconj_free -> ratom_free_members.
 *
 * The concrete solving engines (Phases 1+) are wired on top of this layer; this
 * header is the shared vocabulary between reduce.c, reduce_atom.c and the
 * per-domain solvers.
 */
#ifndef REDUCE_FORM_H
#define REDUCE_FORM_H

#include "expr.h"
#include <stdbool.h>

/* Relation kind of a canonical atom.  Greater/GreaterEqual never appear -- the
 * canonicaliser swaps sides to express them as R_LT / R_LE. */
typedef enum { R_EQ, R_NE, R_LT, R_LE, R_ELEM } RRel;

/* One atom.  For R_EQ/R_NE/R_LT/R_LE it means  poly REL 0.
 * For R_ELEM it means  Element[poly, elem_dom]. */
typedef struct {
    Expr*  poly;       /* owned; canonical LHS-RHS with RHS forced to 0     */
    RRel   rel;
    Expr*  elem_dom;   /* owned; the domain symbol for R_ELEM, else NULL    */
    /* Optional "solved" display override, owned or NULL.  When non-NULL the
     * atom EMITS as this verbatim relation Expr (e.g. Equal[x, b/a]) instead of
     * the canonical `poly REL 0` form, so a solution reads as `x == value`
     * rather than `x - value == 0`.  `poly`/`rel` still carry the semantic
     * content used for decision/dedup, so `display` never affects reasoning. */
    Expr*  display;
    /* True when canonicalisation (Numerator[Together[lhs-rhs]]) cleared a
     * denominator that is NOT a nonzero constant -- i.e. the relation is a
     * rational function whose polynomial numerator has lost the sign
     * information of the (variable) denominator and its poles.  The real
     * engines (sign diagram, Fourier-Motzkin) must NOT treat such an atom as
     * polynomial; they decline instead, keeping Reduce sound. */
    bool   nonconst_denom;
    /* Cheap cached classification filled by reduce_atom_canonicalize.  The
     * degree / linearity fields are Phase-2 routing hints; Phase 0 leaves them
     * at the sentinels below. */
    int    main_var;   /* index into the ambient var array, or -1           */
    int    deg_main;   /* degree in main_var, or -1 when not computed        */
    bool   is_linear;  /* total degree <= 1 in all reduce vars              */
} RAtom;

/* A conjunction (logical And) of atoms. */
typedef struct {
    RAtom* a;
    int    n, cap;
    bool   is_false;   /* a proven contradiction: the whole conjunction is False */
} RConj;

/* A disjunction (logical Or) of conjunctions -- disjunctive normal form. */
typedef struct {
    RConj** c;
    int     n, cap;
    bool    is_true;   /* a proven tautology: the whole form is True */
} RForm;

/* ---- RAtom ------------------------------------------------------------- */

/* Free the owned members of `a` (poly, elem_dom).  Does not free `a` itself
 * (RAtoms live by value inside an RConj's array). */
void  ratom_free_members(RAtom* a);
/* Deep copy (owned poly / elem_dom duplicated). */
RAtom ratom_dup(const RAtom* a);
/* Structural equality of the relation and polynomial (metadata ignored). */
bool  ratom_eq(const RAtom* x, const RAtom* y);
/* Logical negation of a single atom, as an atom (De Morgan leaf):
 *   !(p==0) -> p!=0,  !(p!=0) -> p==0,
 *   !(p<0)  -> -p<=0, !(p<=0) -> -p<0.
 * R_ELEM negation is unsupported in Phase 0: sets *ok=false. */
RAtom ratom_negate(const RAtom* a, bool* ok);

/* ---- RConj ------------------------------------------------------------- */

RConj* rconj_new(void);
void   rconj_free(RConj* c);
RConj* rconj_dup(const RConj* c);
/* Append `a` (taking ownership of its members).  Skips an atom already present
 * (structural dedup). */
void   rconj_push(RConj* c, RAtom a);

/* ---- RForm ------------------------------------------------------------- */

RForm* rform_new(void);
void   rform_free(RForm* f);
RForm* rform_true(void);      /* the tautology  */
RForm* rform_false(void);     /* the contradiction */
RForm* rform_single_atom(RAtom a);   /* one conjunction with one atom */
/* Append a conjunction (taking ownership of `c`). */
void   rform_add_conj(RForm* f, RConj* c);

/* Logical Or / And of two forms.  Both operands are CONSUMED; a freshly-owned
 * form is returned. */
RForm* rform_or(RForm* a, RForm* b);
RForm* rform_and(RForm* a, RForm* b);

/* In-place cleanup: decide constant atoms (via evaluate), drop `is_false`
 * conjunctions, promote an empty conjunction to `is_true`, and dedup identical
 * conjunctions.  `vars`/`nv` are the ambient reduce variables (unused in the
 * Phase-0 constant pass but threaded for later bound-merging). */
void   rform_simplify(RForm* f, Expr** vars, int nv);

/* Emit the form as a Mathematica-form Expr (And/Or/relational atoms), then a
 * single evaluate() pass to flatten.  True -> Symbol True, False -> Symbol
 * False.  Returns a freshly-owned Expr. */
Expr*  rform_to_expr(const RForm* f, Expr** vars, int nv);

/* ---- Atom canonicalisation / form construction (reduce_atom.c) --------- */

/* Canonicalise a single relational Expr (Equal/Unequal/Less/LessEqual/Greater/
 * GreaterEqual/Element) into an RAtom.  Sets *ok=false for an unsupported head.
 * `vars`/`nv` are used only to fill the `main_var` hint. */
RAtom  reduce_atom_canonicalize(const Expr* rel, Expr** vars, int nv, bool* ok);

/* Decide a *constant* atom by evaluating its relation: returns 1 (true), 0
 * (false), or -1 (undecided -- symbolic, parametric, or genuinely unknown). */
int    reduce_atom_decide_constant(const RAtom* a);

/* Build a "solved" equation atom  var == value:  canonical poly = var - value
 * (numerator form), rel R_EQ, and a `display` of Equal[var, value] so it emits
 * in solved form.  `var` and `value` are BORROWED (copied internally). */
RAtom  ratom_solved(const Expr* var, const Expr* value, Expr** vars, int nv);

/* Parse a logical combination (&&, ||, !, Implies, binary Xor, chained
 * Inequality, and relational leaves) into a DNF RForm.  Sets *ok=false when the
 * input contains a construct this layer cannot normalise yet; the returned form
 * is then meaningless and must be freed by the caller.  Returns a freshly-owned
 * form. */
RForm* reduce_form_from_expr(const Expr* e, Expr** vars, int nv, bool* ok);

#endif /* REDUCE_FORM_H */
