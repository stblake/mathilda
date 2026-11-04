/* test_match_order.c
 *
 * Enumeration order of the pattern matcher (ReplaceList, and so which match
 * Replace / a rule takes first) for Orderless, Flat + Orderless and ordered
 * heads, against the reference kernel.  The rule, documented at "Enumeration
 * order" in src/match.c, for the matches of one level (one call's arguments):
 *   1. outermost, the Optional defaults of the level's slots: all matches
 *      with no default first, then those with the last Optional slot
 *      defaulted, ..., the earliest Optional slot counting most;
 *   2. then which elements go to which slot, slots in canonical pattern
 *      order with the leftover-takers (sequences, Flat blanks, and an
 *      Optional wrapping one) last -- for an ordered head, the sequence split;
 *   3. innermost, for that assignment, each slot's own internal alternatives
 *      (orientations of a nested Orderless call, nested splits, Alternatives,
 *      an inner Optional), slot by slot in slot order.
 * Expected values computed by the reference kernel (wolframscript).
 * Not covered here, deliberately: the reference wraps a single operand of a
 * Flat head as fo[a] in some matches (ReplaceList[fo[a, b, c], fo[x_, y_]]),
 * and evaluates a /; inside one slot of an Orderless product without the
 * other slots' bindings; Mathilda differs on both.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

/* lhs and rhs evaluate to the same expression */
static void same(const char* lhs, const char* rhs) {
    char buf[8192];
    snprintf(buf, sizeof buf, "(%s) === (%s)", lhs, rhs);
    assert_eval_eq(buf, "True", 0);
}

static void setup(void) {
    same("SetAttributes[o, Orderless]; SetAttributes[fo, {Flat, Orderless}]; "
         "SetAttributes[Aa, Orderless]; SetAttributes[Zz, Orderless]; "
         "SetAttributes[Bb, Orderless]; SetAttributes[fl, Flat]; True", "True");
}

/* The original battery: Orderless sequences and blanks, Optional leftover
 * (r_. s[a_], r_. s[a_] s[b_]) and the loop nesting of a nested Orderless
 * slot (r_ ch[a_, ii_] Aa[ii_, jj_]). */
static void test_orderless_battery(void) {
    same("ReplaceList[o[a, b, c], o[x_, y__] :> {x, {y}}]",
         "{{a, {b, c}}, {b, {a, c}}, {c, {a, b}}}");
    same("ReplaceList[o[a, b, c], o[x__, y_] :> {{x}, y}]",
         "{{{b, c}, a}, {{a, c}, b}, {{a, b}, c}}");
    same("ReplaceList[o[a, b, c], o[x_, y_, z_] :> {x, y, z}]",
         "{{a, b, c}, {a, c, b}, {b, a, c}, {b, c, a}, {c, a, b}, {c, b, a}}");
    same("ReplaceList[o[a, b, c], o[x___, y_] :> {{x}, y}]",
         "{{{b, c}, a}, {{a, c}, b}, {{a, b}, c}}");
    same("ReplaceList[o[a, b, c], o[x_, y___] :> {x, {y}}]",
         "{{a, {b, c}}, {b, {a, c}}, {c, {a, b}}}");
    same("ReplaceList[o[a, b, c, d], o[x__, y__] :> {{x}, {y}}]",
         "{{{a}, {b, c, d}}, {{b}, {a, c, d}}, {{c}, {a, b, d}}, {{d}, {a, b, c}}, {{a, b}, {c, d}}, {{a, c}, {b, d}}, {{a, d}, {b, c}}, {{b, c}, {a, d}}, {{b, d}, {a, c}}, {{c, d}, {a, b}}, {{a, b, c}, {d}}, {{a, b, d}, {c}}, {{a, c, d}, {b}}, {{b, c, d}, {a}}}");
    same("ReplaceList[o[s[1], s[2], s[3]], o[s[a_], r___] :> {a, {r}}]",
         "{{1, {s[2], s[3]}}, {2, {s[1], s[3]}}, {3, {s[1], s[2]}}}");
    same("ReplaceList[o[s[1], s[2], s[3]], o[r___, s[a_]] :> {a, {r}}]",
         "{{1, {s[2], s[3]}}, {2, {s[1], s[3]}}, {3, {s[1], s[2]}}}");
    same("ReplaceList[o[s[1], s[2], s[3]], o[r_, s[a_], s[b_]] :> {r, a, b}]",
         "{{s[1], 2, 3}, {s[1], 3, 2}, {s[2], 1, 3}, {s[2], 3, 1}, {s[3], 1, 2}, {s[3], 2, 1}}");
    same("ReplaceList[s[1] s[2] s[3], r_ s[a_] :> {r, a}]",
         "{{s[2]*s[3], 1}, {s[1]*s[3], 2}, {s[1]*s[2], 3}}");
    same("ReplaceList[s[1] s[2] s[3], r_. s[a_] :> {r, a}]",
         "{{s[2]*s[3], 1}, {s[1]*s[3], 2}, {s[1]*s[2], 3}}");
    same("ReplaceList[s[1] s[2] s[3] s[4], r_. s[a_] s[b_] :> {r, a, b}]",
         "{{s[3]*s[4], 1, 2}, {s[2]*s[4], 1, 3}, {s[2]*s[3], 1, 4}, {s[3]*s[4], 2, 1}, {s[1]*s[4], 2, 3}, {s[1]*s[3], 2, 4}, {s[2]*s[4], 3, 1}, {s[1]*s[4], 3, 2}, {s[1]*s[2], 3, 4}, {s[2]*s[3], 4, 1}, {s[1]*s[3], 4, 2}, {s[1]*s[2], 4, 3}}");
    same("ReplaceList[s[1] s[2] s[3], s[a_] s[b_] :> {a, b}]",
         "{}");
    same("ReplaceList[s[1] s[2] s[3], r__ s[a_] :> {{r}, a}]",
         "{{{s[2], s[3]}, 1}, {{s[1], s[3]}, 2}, {{s[1], s[2]}, 3}}");
    same("ReplaceList[o[a, b, c], o[x_, y_] :> {x, y}]",
         "{}");
    same("ReplaceList[Aa[1, 2] g[1] g[2], Aa[i_, j_] g[i_] r_ :> {i, j, r}]",
         "{{1, 2, g[2]}, {2, 1, g[1]}}");
    same("ReplaceList[ch[e, j] ch[x, i] Aa[i, j], r_ ch[a_, ii_] Aa[ii_, jj_] :> {ii, r}]",
         "{{j, ch[x, i]}, {i, ch[e, j]}}");
    same("ReplaceList[ch[e, j] ch[x, i] Aa[i, j], ch[a_, ii_] Aa[ii_, jj_] r_ :> {ii, r}]",
         "{{j, ch[x, i]}, {i, ch[e, j]}}");
    same("ReplaceList[Aa[L[4], L[5]] hd[1, L[5]] hd[3, L[4]], Aa[L[z_], L[x_]] f_[u__] :> {z, x, f[u]}]",
         "{}");
    same("ReplaceList[Aa[L[4], L[5]] hd[1, L[5]] hd[3, L[4]], Aa[L[z_], L[x_]] f_[u__] /; !FreeQ[f[u], L[z]] :> {z, x, f[u]}]",
         "{}");
    same("ReplaceList[o[1, 2, 3], o[x_?OddQ, y__] :> {x, {y}}]",
         "{{1, {2, 3}}, {3, {1, 2}}}");
    same("ReplaceList[o[a, a, b], o[x_, x_, y_] :> {x, y}]",
         "{{a, b}, {a, b}}");
    same("ReplaceList[s[1] s[2] s[3] s[4], s[a_] s[b_] r__ :> {a, b, {r}}]",
         "{{1, 2, {s[3], s[4]}}, {1, 3, {s[2], s[4]}}, {1, 4, {s[2], s[3]}}, {2, 1, {s[3], s[4]}}, {2, 3, {s[1], s[4]}}, {2, 4, {s[1], s[3]}}, {3, 1, {s[2], s[4]}}, {3, 2, {s[1], s[4]}}, {3, 4, {s[1], s[2]}}, {4, 1, {s[2], s[3]}}, {4, 2, {s[1], s[3]}}, {4, 3, {s[1], s[2]}}}");
    same("ReplaceList[Plus[a, b, c, d], x_ + y_ :> {x, y}]",
         "{{a, b + c + d}, {b, a + c + d}, {c, a + b + d}, {d, a + b + c}, {a + b, c + d}, {a + c, b + d}, {a + d, b + c}, {b + c, a + d}, {b + d, a + c}, {c + d, a + b}, {a + b + c, d}, {a + b + d, c}, {a + c + d, b}, {b + c + d, a}}");
    same("ReplaceList[Plus[a, b, c], x_. + y_ :> {x, y}]",
         "{{a, b + c}, {b, a + c}, {c, a + b}, {a + b, c}, {a + c, b}, {b + c, a}, {0, a + b + c}}");
}

/* Optional defaults (outermost), nested alternatives in several slots,
 * Alternatives, nested splits, inner Optional, conditions, ordered heads. */
static void test_model_probes(void) {
    same("Replace[s[1] s[2], r_. s[a_] s[b_] :> {r, a, b}]",
         "{1, 1, 2}");
    same("ReplaceList[s[1] s[2], r_. s[a_] s[b_] :> {r, a, b}]",
         "{{1, 1, 2}, {1, 2, 1}}");
    same("ReplaceList[s[1] s[2] s[3], r_. s[a_] s[b_] :> {r, a, b}]",
         "{{s[3], 1, 2}, {s[2], 1, 3}, {s[3], 2, 1}, {s[1], 2, 3}, {s[2], 3, 1}, {s[1], 3, 2}}");
    same("ReplaceList[o[s[1], s[2]], o[r_:0, s[a_]] :> {r, a}]",
         "{{s[1], 2}, {s[2], 1}}");
    same("ReplaceList[o[s[1], s[2]], o[r_:0, s[a_], s[b_]] :> {r, a, b}]",
         "{{0, 1, 2}, {0, 2, 1}}");
    same("ReplaceList[o[s[1], s[2], s[3]], o[r_:0, s[a_], s[b_]] :> {r, a, b}]",
         "{{s[1], 2, 3}, {s[1], 3, 2}, {s[2], 1, 3}, {s[2], 3, 1}, {s[3], 1, 2}, {s[3], 2, 1}}");
    same("ReplaceList[s[1] + s[2] + s[3], r_. + s[a_] :> {r, a}]",
         "{{s[2] + s[3], 1}, {s[1] + s[3], 2}, {s[1] + s[2], 3}}");
    same("ReplaceList[s[1] s[2] s[3], r_. q_. s[a_] :> {r, q, a}]",
         "{{s[3], s[2], 1}, {s[2], s[3], 1}, {s[3], s[1], 2}, {s[1], s[3], 2}, {s[2], s[1], 3}, {s[1], s[2], 3}, {1, s[2]*s[3], 1}, {1, s[1]*s[3], 2}, {1, s[1]*s[2], 3}, {s[2]*s[3], 1, 1}, {s[1]*s[3], 1, 2}, {s[1]*s[2], 1, 3}}");
    same("ReplaceList[s[1] s[2] s[3], r_ q_. s[a_] :> {r, q, a}]",
         "{{s[3], s[2], 1}, {s[2], s[3], 1}, {s[3], s[1], 2}, {s[1], s[3], 2}, {s[2], s[1], 3}, {s[1], s[2], 3}, {s[2]*s[3], 1, 1}, {s[1]*s[3], 1, 2}, {s[1]*s[2], 1, 3}}");
    same("ReplaceList[Aa[1, 2] g[5] g[6] Zz[3, 4], Aa[i_, j_] g[x_] Zz[k_, l_] r_ :> {i, j, x, k, l, r}]",
         "{{1, 2, 5, 3, 4, g[6]}, {1, 2, 5, 4, 3, g[6]}, {2, 1, 5, 3, 4, g[6]}, {2, 1, 5, 4, 3, g[6]}, {1, 2, 6, 3, 4, g[5]}, {1, 2, 6, 4, 3, g[5]}, {2, 1, 6, 3, 4, g[5]}, {2, 1, 6, 4, 3, g[5]}}");
    same("ReplaceList[Aa[1, 2] g[5] Zz[3, 4], Aa[i_, j_] g[x_] Zz[k_, l_] :> {i, j, x, k, l}]",
         "{{1, 2, 5, 3, 4}, {1, 2, 5, 4, 3}, {2, 1, 5, 3, 4}, {2, 1, 5, 4, 3}}");
    same("ReplaceList[Aa[1, 2] Bb[5, 6] Zz[3, 4], Aa[i_, j_] Bb[x_, y_] Zz[k_, l_] :> {i, j, x, y, k, l}]",
         "{{1, 2, 5, 6, 3, 4}, {1, 2, 5, 6, 4, 3}, {1, 2, 6, 5, 3, 4}, {1, 2, 6, 5, 4, 3}, {2, 1, 5, 6, 3, 4}, {2, 1, 5, 6, 4, 3}, {2, 1, 6, 5, 3, 4}, {2, 1, 6, 5, 4, 3}}");
    same("ReplaceList[p[1, 2] q[5] q[6], (p[x_, _] | p[_, x_]) q[y_] r_ :> {x, y, r}]",
         "{{1, 5, q[6]}, {2, 5, q[6]}, {1, 6, q[5]}, {2, 6, q[5]}}");
    same("ReplaceList[p[1, 2] q[5] q[6], q[y_] (p[x_, _] | p[_, x_]) r_ :> {x, y, r}]",
         "{{1, 5, q[6]}, {2, 5, q[6]}, {1, 6, q[5]}, {2, 6, q[5]}}");
    same("ReplaceList[f[1, 2, 3] q[5] q[6], f[x__, y__] q[z_] r_ :> {{x}, {y}, z, r}]",
         "{{{1}, {2, 3}, 5, q[6]}, {{1, 2}, {3}, 5, q[6]}, {{1}, {2, 3}, 6, q[5]}, {{1, 2}, {3}, 6, q[5]}}");
    same("ReplaceList[g[1] q[5] q[6], g[x_:0, y_:0] q[z_] r_ :> {x, y, z, r}]",
         "{{1, 0, 5, q[6]}, {0, 1, 5, q[6]}, {1, 0, 6, q[5]}, {0, 1, 6, q[5]}}");
    same("ReplaceList[Aa[1, 2] g[1] g[2] g[3], Aa[i_, j_] g[k_] r_ /; k > i :> {i, j, k, r}]",
         "{{1, 2, 2, g[1]*g[3]}, {1, 2, 3, g[1]*g[2]}, {2, 1, 3, g[1]*g[2]}}");
    same("ReplaceList[Aa[1, 2] g[1] g[2] g[3], Aa[i_, j_] g[k_ /; k > 1] r_ :> {i, j, k, r}]",
         "{{1, 2, 2, g[1]*g[3]}, {2, 1, 2, g[1]*g[3]}, {1, 2, 3, g[1]*g[2]}, {2, 1, 3, g[1]*g[2]}}");
    same("ReplaceList[Aa[1, 2] g[1] g[2] g[3], Aa[i_, j_] g[k_] r_ /; k > j :> {i, j, k, r}]",
         "{{2, 1, 2, g[1]*g[3]}, {1, 2, 3, g[1]*g[2]}, {2, 1, 3, g[1]*g[2]}}");
    same("ReplaceList[Aa[1, 2] Bb[2, 3] Bb[1, 4], Aa[i_, j_] Bb[j_, k_] r_ :> {i, j, k, r}]",
         "{{2, 1, 4, Bb[2, 3]}, {1, 2, 3, Bb[1, 4]}}");
    same("ReplaceList[f[Aa[1, 2], 3, 4, 5], f[Aa[i_, j_], x__, y__] :> {i, {x}, {y}}]",
         "{{1, {3}, {4, 5}}, {2, {3}, {4, 5}}, {1, {3, 4}, {5}}, {2, {3, 4}, {5}}}");
    same("ReplaceList[f[Aa[1, 2], Aa[3, 4]], f[Aa[i_, j_], Aa[k_, l_]] :> {i, j, k, l}]",
         "{{1, 2, 3, 4}, {1, 2, 4, 3}, {2, 1, 3, 4}, {2, 1, 4, 3}}");
    same("ReplaceList[o[Aa[1, 2], g[5], g[6]], o[Aa[i_, j_], g[x_], r_] :> {i, x, r}]",
         "{{1, 5, g[6]}, {2, 5, g[6]}, {1, 6, g[5]}, {2, 6, g[5]}}");
    same("ReplaceList[o[Aa[1, 2], g[5], g[6]], o[Aa[i_, j_], g[x_], r___] :> {i, x, {r}}]",
         "{{1, 5, {g[6]}}, {2, 5, {g[6]}}, {1, 6, {g[5]}}, {2, 6, {g[5]}}}");
    same("ReplaceList[Aa[1, 2] Aa[3, 4] g[5], Aa[i_, j_] Aa[k_, l_] r_ :> {i, j, k, l, r}]",
         "{{1, 2, 3, 4, g[5]}, {1, 2, 4, 3, g[5]}, {2, 1, 3, 4, g[5]}, {2, 1, 4, 3, g[5]}, {3, 4, 1, 2, g[5]}, {3, 4, 2, 1, g[5]}, {4, 3, 1, 2, g[5]}, {4, 3, 2, 1, g[5]}}");
    same("ReplaceList[Aa[1, 2] Aa[3, 4] g[5] g[6], Aa[i_, j_] g[x_] r_ :> {i, j, x, r}]",
         "{{1, 2, 5, Aa[3, 4]*g[6]}, {2, 1, 5, Aa[3, 4]*g[6]}, {1, 2, 6, Aa[3, 4]*g[5]}, {2, 1, 6, Aa[3, 4]*g[5]}, {3, 4, 5, Aa[1, 2]*g[6]}, {4, 3, 5, Aa[1, 2]*g[6]}, {3, 4, 6, Aa[1, 2]*g[5]}, {4, 3, 6, Aa[1, 2]*g[5]}}");
    same("ReplaceList[Aa[1, 2] g[5] g[6] g[7], Aa[i_, j_] g[x_] g[y_] r_ :> {i, x, y, r}]",
         "{{1, 5, 6, g[7]}, {2, 5, 6, g[7]}, {1, 5, 7, g[6]}, {2, 5, 7, g[6]}, {1, 6, 5, g[7]}, {2, 6, 5, g[7]}, {1, 6, 7, g[5]}, {2, 6, 7, g[5]}, {1, 7, 5, g[6]}, {2, 7, 5, g[6]}, {1, 7, 6, g[5]}, {2, 7, 6, g[5]}}");
    same("ReplaceList[w[Aa[1, 2]] g[5] g[6], w[Aa[i_, j_]] g[x_] r_ :> {i, x, r}]",
         "{{1, 5, g[6]}, {2, 5, g[6]}, {1, 6, g[5]}, {2, 6, g[5]}}");
    same("ReplaceList[Aa[1, 2] g[5] g[6], Aa[i_, j_] g[x_] r__ :> {i, x, {r}}]",
         "{{1, 5, {g[6]}}, {2, 5, {g[6]}}, {1, 6, {g[5]}}, {2, 6, {g[5]}}}");
    same("ReplaceList[Aa[1, 2] g[5] g[6], Aa[i_, j_] g[x_] r___ :> {i, x, {r}}]",
         "{{1, 5, {g[6]}}, {2, 5, {g[6]}}, {1, 6, {g[5]}}, {2, 6, {g[5]}}}");
    same("ReplaceList[Aa[1, 2] g[5] g[6], Aa[i_, j_] r_. g[x_] :> {i, x, r}]",
         "{{1, 5, g[6]}, {2, 5, g[6]}, {1, 6, g[5]}, {2, 6, g[5]}}");
    same("ReplaceList[Aa[1, 2] + g[5] + g[6], Aa[i_, j_] + g[x_] + r_ :> {i, x, r}]",
         "{{1, 5, g[6]}, {2, 5, g[6]}, {1, 6, g[5]}, {2, 6, g[5]}}");
}

/* Linked variables across nested Orderless slots, first-match Replace, ordered
 * and nested levels, several Optionals, MatchQ. */
static void test_model_probes_2(void) {
    same("ReplaceList[Aa[1, 2] Aa[2, 3] Aa[3, 1] g[1], Aa[i_, j_] Aa[j_, k_] r_ :> {i, j, k, r}]",
         "{{2, 1, 3, Aa[2, 3]*g[1]}, {1, 2, 3, Aa[1, 3]*g[1]}, {3, 1, 2, Aa[2, 3]*g[1]}, {1, 3, 2, Aa[1, 2]*g[1]}, {3, 2, 1, Aa[1, 3]*g[1]}, {2, 3, 1, Aa[1, 2]*g[1]}}");
    same("ReplaceList[Aa[1, 2] Aa[2, 3] Aa[3, 1], Aa[i_, j_] Aa[j_, k_] Aa[k_, l_] :> {i, j, k, l}]",
         "{{2, 1, 3, 2}, {1, 2, 3, 1}, {3, 1, 2, 3}, {1, 3, 2, 1}, {3, 2, 1, 3}, {2, 3, 1, 2}}");
    same("ReplaceList[Aa[1, 2] g[2] g[1] h[3], Aa[i_, j_] g[j_] r__ :> {i, j, {r}}]",
         "{{2, 1, {g[2], h[3]}}, {1, 2, {g[1], h[3]}}}");
    same("Replace[Aa[1, 2] g[2] g[1] h[3], Aa[i_, j_] g[j_] r__ :> {i, j, {r}}]",
         "{2, 1, {g[2], h[3]}}");
    same("Replace[ch[e, j] ch[x, i] Aa[i, j], r_ ch[a_, ii_] Aa[ii_, jj_] :> {ii, r}]",
         "{j, ch[x, i]}");
    same("ReplaceList[f[Aa[1, 2], 5], f[Aa[i_, j_], x_:0, y_:0] :> {i, x, y}]",
         "{{1, 5, 0}, {2, 5, 0}, {1, 0, 5}, {2, 0, 5}}");
    same("ReplaceList[f[Aa[1, 2], 5, 6], f[Aa[i_, j_], x___, y___] :> {i, {x}, {y}}]",
         "{{1, {}, {5, 6}}, {2, {}, {5, 6}}, {1, {5}, {6}}, {2, {5}, {6}}, {1, {5, 6}, {}}, {2, {5, 6}, {}}}");
    same("ReplaceList[o[Aa[1, 2], Zz[3, 4], 5, 6], o[Aa[i_, j_], Zz[k_, l_], x_, y_] :> {i, k, x, y}]",
         "{{1, 3, 5, 6}, {1, 4, 5, 6}, {2, 3, 5, 6}, {2, 4, 5, 6}, {1, 3, 6, 5}, {1, 4, 6, 5}, {2, 3, 6, 5}, {2, 4, 6, 5}}");
    same("ReplaceList[o[Aa[1, 2], Zz[3, 4], 5], o[Zz[k_, l_], Aa[i_, j_], x_] :> {i, k, x}]",
         "{{1, 3, 5}, {1, 4, 5}, {2, 3, 5}, {2, 4, 5}}");
    same("ReplaceList[u[Aa[1, 2] g[5] g[6]], u[Aa[i_, j_] g[x_] r_] :> {i, x, r}]",
         "{{1, 5, g[6]}, {2, 5, g[6]}, {1, 6, g[5]}, {2, 6, g[5]}}");
    same("ReplaceList[{Aa[1, 2] g[5] g[6], 7}, {Aa[i_, j_] g[x_] r_, y_} :> {i, x, r, y}]",
         "{{1, 5, g[6], 7}, {2, 5, g[6], 7}, {1, 6, g[5], 7}, {2, 6, g[5], 7}}");
    same("ReplaceList[Aa[1, 2] g[5] g[6], Aa[i_, j_] g[x_] r_ /; x > 5 :> {i, x, r}]",
         "{{1, 6, g[5]}, {2, 6, g[5]}}");
    same("ReplaceList[Aa[1, 2] g[5] g[6], (Aa[i_, j_] g[x_] r_) /; i > 1 :> {i, x, r}]",
         "{{2, 5, g[6]}, {2, 6, g[5]}}");
    same("ReplaceList[Aa[1, 2] Aa[1, 3] g[5], Aa[1, j_] g[x_] r_ :> {j, x, r}]",
         "{{2, 5, Aa[1, 3]}, {3, 5, Aa[1, 2]}}");
    same("ReplaceList[Aa[x, y] + Aa[y, z] + c, Aa[a_, b_] + Aa[b_, d_] + r_ :> {a, b, d, r}]",
         "{{x, y, z, c}, {z, y, x, c}}");
    same("ReplaceList[s[1] s[2] s[3] s[4], r_. q_. s[a_] s[b_] :> {r, q, a, b}][[;; 6]]",
         "{{s[4], s[3], 1, 2}, {s[3], s[4], 1, 2}, {s[4], s[2], 1, 3}, {s[2], s[4], 1, 3}, {s[3], s[2], 1, 4}, {s[2], s[3], 1, 4}}");
    same("ReplaceList[s[1] s[2] s[3], r__ q_. s[a_] :> {{r}, q, a}]",
         "{{{s[3]}, s[2], 1}, {{s[2]}, s[3], 1}, {{s[3]}, s[1], 2}, {{s[1]}, s[3], 2}, {{s[2]}, s[1], 3}, {{s[1]}, s[2], 3}, {{s[2], s[3]}, 1, 1}, {{s[1], s[3]}, 1, 2}, {{s[1], s[2]}, 1, 3}}");
    same("ReplaceList[s[1] + s[2] + s[3], r_. + q_. + s[a_] :> {r, q, a}]",
         "{{s[3], s[2], 1}, {s[2], s[3], 1}, {s[3], s[1], 2}, {s[1], s[3], 2}, {s[2], s[1], 3}, {s[1], s[2], 3}, {0, s[2] + s[3], 1}, {0, s[1] + s[3], 2}, {0, s[1] + s[2], 3}, {s[2] + s[3], 0, 1}, {s[1] + s[3], 0, 2}, {s[1] + s[2], 0, 3}}");
    same("ReplaceList[a b c, x_. y_. :> {x, y}]",
         "{{a, b*c}, {b, a*c}, {c, a*b}, {a*b, c}, {a*c, b}, {b*c, a}, {a*b*c, 1}, {1, a*b*c}}");
    same("ReplaceList[Aa[1, 2] g[3], Aa[i_, j_] g[x_] r_. :> {i, j, x, r}]",
         "{{1, 2, 3, 1}, {2, 1, 3, 1}}");
    same("ReplaceList[Aa[1, 2] g[3] g[4], Aa[i_, j_] g[x_] r_. :> {i, j, x, r}]",
         "{{1, 2, 3, g[4]}, {2, 1, 3, g[4]}, {1, 2, 4, g[3]}, {2, 1, 4, g[3]}}");
    same("ReplaceList[Aa[1, 2] g[3] g[4], r_. Aa[i_, j_] g[x_] :> {i, j, x, r}]",
         "{{1, 2, 3, g[4]}, {2, 1, 3, g[4]}, {1, 2, 4, g[3]}, {2, 1, 4, g[3]}}");
    same("MatchQ[Aa[1, 2] Aa[2, 3] g[9], Aa[i_, j_] Aa[j_, k_] g[_]]",
         "True");
    same("MatchQ[Aa[1, 2] Aa[4, 3] g[9], Aa[i_, j_] Aa[j_, k_] g[_]]",
         "False");
    same("ReplaceList[w[Aa[1, 2], Aa[3, 4]] g[5] g[6], w[Aa[i_, j_], Aa[k_, l_]] g[x_] r_ :> {i, k, x}]",
         "{{1, 3, 5}, {1, 4, 5}, {2, 3, 5}, {2, 4, 5}, {1, 3, 6}, {1, 4, 6}, {2, 3, 6}, {2, 4, 6}}");
    same("ReplaceList[Aa[p[1], p[2]] q[1] q[2] q[3], Aa[p[i_], p[j_]] q[i_] r__ :> {i, j, {r}}]",
         "{{1, 2, {q[2], q[3]}}, {2, 1, {q[1], q[3]}}}");
}

/* The same order through DownValues and UpValues: the FeynCalc PairContract
 * shape -- the f slot takes hd[1, L[i4]] first and pc's second orientation
 * satisfies the condition before f moves on. */
static void test_definitions(void) {
    same("SetAttributes[pc, Orderless]; "
         "pc /: pc[L[z_], L[x_]] f_[a__] := (f[a] /. L[z] -> L[x]) /; !FreeQ[f[a], L[z]]; "
         "pc[L[i4], L[i5]] hd[1, L[i5]] hd[3, L[i4]]",
         "hd[1, L[i4]]*hd[3, L[i4]]");
    same("pc[L[i4], L[i5]] hd[1, L[i4]] hd[3, L[i5]]", "hd[1, L[i5]]*hd[3, L[i5]]");
    same("pc[L[i5], L[i4]] hd[3, L[i5]] hd[1, L[i4]] hd[0, L[i5]]",
         "hd[0, L[i4]]*hd[1, L[i4]]*hd[3, L[i5]]");
    same("f1[Aa[i_, j_] g[x_] r_] := {i, x, r}; f1[Aa[1, 2] g[5] g[6]]", "{1, 5, g[6]}");
    same("Aa /: Aa[i_, j_] g[j_] r_ := hit[i, j, r]; Aa[1, 2] g[1] g[2] h[3]",
         "hit[2, 1, g[2]*h[3]]");
}

int main(void) {
    symtab_init();
    core_init();
    setup();
    TEST(test_orderless_battery);
    TEST(test_model_probes);
    TEST(test_model_probes_2);
    TEST(test_definitions);
    printf("All match-order tests passed.\n");
    return 0;
}
