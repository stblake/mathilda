/* gb_f4.h -- Faugere's F4 over GF(p) (gb_f4.c).
 *
 * The modular Groebner engine behind gb_flint.c: every GF(p) basis the FLINT
 * engine computes -- Modulus -> p, and the per-prime bases of the modular
 * routes over Q -- comes from here.  Inputs and outputs are nmod_mpoly in one
 * context (ORD_LEX, ORD_DEGLEX or ORD_DEGREVLEX, any word-size prime). */
#ifndef GB_F4_H
#define GB_F4_H

#ifdef USE_FLINT
#include <flint/flint.h>
#include <flint/nmod_mpoly.h>

/* A Groebner basis of <F[0..n)> (F read, never taken): *G (allocated, ng
 * initialised nmod_mpoly, caller clears and flint_frees) is the REDUCED basis,
 * monic, ascending by leading monomial ({1} for the unit ideal).  Returns 1, or 0 when the
 * deadline passed (*expired set; nothing is held) -- or when F4 declines
 * (exponents past 2^20: *expired clear), so the caller's own engine runs. */
int gb_f4_nmod(const nmod_mpoly_struct* F, slong n, const nmod_mpoly_ctx_t ctx,
               nmod_mpoly_struct** G, slong* ng, int* expired);

/* MATHILDA_GB_F4=off, or gb_f4_enable(0), makes gb_f4_nmod decline (A/B
 * tests against the Buchberger engine).  Returns the previous setting. */
int gb_f4_enable(int on);
#endif

#endif /* GB_F4_H */
