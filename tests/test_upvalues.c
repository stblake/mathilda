/* Tests for the UpValues subsystem (trial, gated on UP_VALUES):
 * UpSet (^=), UpSetDelayed (^:=), TagSet (f/:lhs=rhs), TagSetDelayed (f/:lhs:=rhs),
 * TagUnset (f/:lhs=.), UpValues[] reader, the {Up,Down,Own}Values[f]=list setters,
 * Definition[], ?name display, Hold/HoldComplete interaction, scoping, and the
 * efficiency gate counter. Every case mirrors an example from the Wolfram
 * reference the feature was built against.
 *
 * Each test uses a distinct symbol namespace because the symbol table is global
 * state that persists across the sequential TEST() calls. */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

/* Evaluate and discard -- used to install definitions before asserting. */
static void ev(const char* s) {
    Expr* p = parse_expression(s);
    assert(p != NULL);
    Expr* r = evaluate(p);
    expr_free(p);
    expr_free(r);
}

/* ---- UpSet (^=) ---- */
static void test_upset_basic(void) {
    ev("area1[square1] ^= s^2");
    assert_eval_eq("UpValues[square1]", "{HoldPattern[area1[square1]] :> s^2}", 0);
    assert_eval_eq("area1[square1]", "s^2", 0);
    /* not keyed on the outer head */
    assert_eval_eq("UpValues[area1]", "{}", 0);
}

static void test_upset_multi_symbol(void) {
    /* installs on EVERY distinct level-one symbol: a2 (arg) and b2 (arg head) */
    ev("prop2[a2, b2[c2]] ^= value2");
    assert_eval_eq("UpValues[a2]", "{HoldPattern[prop2[a2, b2[c2]]] :> value2}", 0);
    assert_eval_eq("UpValues[b2]", "{HoldPattern[prop2[a2, b2[c2]]] :> value2}", 0);
}

static void test_upset_immediate_rhs(void) {
    /* ^= evaluates its RHS once at definition time */
    ev("rate3[chf3] ^= 2 + 3");
    assert_eval_eq("UpValues[chf3]", "{HoldPattern[rate3[chf3]] :> 5}", 0);
}

/* ---- UpSetDelayed (^:=) ---- */
static void test_upsetdelayed(void) {
    ev("f4[g4[x_]] ^:= h4[x]");
    assert_eval_eq("f4[g4[2]]", "h4[2]", 0);
    assert_eval_eq("f4[hh4[2]]", "f4[hh4[2]]", 0);  /* unrelated head unchanged */
}

static void test_upsetdelayed_arg_head(void) {
    ev("area5[sq5[s_]] ^:= s^2");
    assert_eval_eq("area5[sq5[3]]", "9", 0);
}

/* ---- TagSet / TagSetDelayed (f /: ...) ---- */
static void test_tagsetdelayed_upvalue(void) {
    ev("g6 /: f6[g6[x_]] := h6[x]");
    assert_eval_eq("f6[g6[2]]", "h6[2]", 0);
    assert_eval_eq("UpValues[g6]", "{HoldPattern[f6[g6[x_]]] :> h6[x]}", 0);
}

static void test_tagset_downvalue_redundant(void) {
    /* tag is the head of lhs -> DownValue (tag redundant) */
    ev("f7 /: f7[x_] := x^2");
    assert_eval_eq("DownValues[f7]", "{HoldPattern[f7[x_]] :> x^2}", 0);
    assert_eval_eq("UpValues[f7]", "{}", 0);
    assert_eval_eq("f7[4]", "16", 0);
}

static void test_tagset_ownvalue_redundant(void) {
    /* tag is lhs itself -> OwnValue (tag redundant) */
    ev("x8 /: x8 = 7");
    assert_eval_eq("OwnValues[x8]", "{HoldPattern[x8] :> 7}", 0);
    assert_eval_eq("x8", "7", 0);
}

static void test_tagset_blank_head(void) {
    /* tag appears as the head of a Blank: a_mod -> installs on the blank's head */
    ev("mod9 /: a_mod9 + b_mod9 := modPlus9[a, b]");
    assert_eval_eq("UpValues[mod9]", "{HoldPattern[a_mod9 + b_mod9] :> modPlus9[a, b]}", 0);
    assert_eval_eq("mod9[1] + mod9[2]", "modPlus9[mod9[1], mod9[2]]", 0);
}

static void test_tagset_not_found(void) {
    /* tag not in lhs: message (suppressed) + yields rhs, installs nothing */
    assert_eval_eq("Quiet[zz10 /: aa10[bb10] = 1]", "1", 0);
    assert_eval_eq("UpValues[zz10]", "{}", 0);
}

/* ---- modular arithmetic end-to-end (the showcase) ---- */
static void test_modular_arithmetic(void) {
    ev("mod11 /: mod11[a_, p_] + mod11[b_, p_] := mod11[Mod[a + b, p], p]");
    ev("mod11 /: i_Integer mod11[a_, p_] := mod11[Mod[i a, p], p]");
    assert_eval_eq("mod11[2, 5] + 3 mod11[3, 5] - mod11[1, 5]", "mod11[0, 5]", 0);
}

/* ---- TagUnset (f /: lhs =.) ---- */
static void test_tagunset(void) {
    ev("h12 /: f12[h12[x_]] := fh12[x]");
    assert_eval_eq("f12[h12[5]]", "fh12[5]", 0);
    ev("h12 /: f12[h12[x_]] =.");
    assert_eval_eq("f12[h12[5]]", "f12[h12[5]]", 0);
    assert_eval_eq("UpValues[h12]", "{}", 0);
}

static void test_plain_unset_keeps_upvalue(void) {
    /* plain `lhs =.` does NOT remove an upvalue (the tag is required) */
    ev("h13 /: f13[h13[x_]] := fh13[x]");
    ev("f13[h13[x_]] =.");
    assert_eval_eq("f13[h13[7]]", "fh13[7]", 0);
}

/* ---- dispatch precedence: upvalue before downvalue ---- */
static void test_upvalue_before_downvalue(void) {
    ev("f14[g14[x_]] ^:= hh14[x]");
    ev("f14[x_] := 1");
    assert_eval_eq("f14[g14[2]]", "hh14[2]", 0);   /* upvalue wins */
    assert_eval_eq("f14[other14]", "1", 0);        /* downvalue for the rest */
}

/* ---- all-heads upvalue + Hold / HoldComplete ---- */
static void test_allheads_hold(void) {
    ev("_[g15[x_]] ^= 1");
    assert_eval_eq("UpValues[g15]", "{HoldPattern[_[g15[x_]]] :> 1}", 0);
    /* fires even under HoldAll (Hold) ... */
    assert_eval_eq("Hold[g15[anything]]", "1", 0);
    /* ... but NOT under HoldAllComplete (HoldComplete) */
    assert_eval_eq("HoldComplete[g15[anything]]", "HoldComplete[g15[anything]]", 0);
}

/* ---- UpValues[] reader ---- */
static void test_upvalues_reader(void) {
    assert_eval_eq("UpValues[neverdef16]", "{}", 0);          /* undefined symbol -> {} */
    ev("g16 /: wrap16[g16] := done16");
    assert_eval_eq("UpValues[\"g16\"]", "{HoldPattern[wrap16[g16]] :> done16}", 0); /* string ok */
    /* string naming a non-existent symbol issues a message (suppressed) */
    assert_eval_eq("Quiet[Head[UpValues[\"neverdef16b\"]]]", "UpValues", 0);
}

/* ---- {Up,Down,Own}Values[f] = list setters ---- */
static void test_upvalues_setter(void) {
    ev("UpValues[h17] = {h17[x_] h17[y_] :> htimes17[x, y]}");
    assert_eval_eq("h17[a] h17[b]", "htimes17[a, b]", 0);
}

static void test_upvalues_copy(void) {
    ev("g18 /: g18[x_] + g18[y_] := plus18[x, y]");
    ev("UpValues[h18] = UpValues[g18] /. g18 -> h18");
    assert_eval_eq("h18[a] + h18[b]", "plus18[a, b]", 0);
}

static void test_upvalues_reorder(void) {
    /* equal-specificity rules: insertion order is the tie-break, so reversing the
     * list changes which fires first */
    ev("x19 /: x19 + y_ /; y > -2 := fpos19[y]");
    ev("x19 /: x19 + y_ /; y < 2 := gpos19[y]");
    assert_eval_eq("x19 + 1", "fpos19[1]", 0);
    ev("UpValues[x19] = Reverse[UpValues[x19]]");
    assert_eval_eq("x19 + 1", "gpos19[1]", 0);
}

static void test_downvalues_setter(void) {
    ev("DownValues[dv20] = {HoldPattern[dv20[x_]] :> x + 100}");
    assert_eval_eq("dv20[5]", "105", 0);
}

static void test_ownvalues_setter(void) {
    ev("OwnValues[ov21] = {HoldPattern[ov21] :> 999}");
    assert_eval_eq("ov21", "999", 0);
}

/* ---- Definition[] ---- */
static void test_definition_inert_fullform(void) {
    ev("g22 /: f22[g22[x_]] := h22[x]");
    /* Definition[g] is inert: its FullForm is the object itself, not its
     * rendered body (is_fullform=1 prints via the FullForm printer). */
    assert_eval_eq("Definition[g22]", "Definition[g22]", 1);
}

/* ---- scoping: Block restores upvalues ---- */
static void test_block_restores_upvalues(void) {
    ev("g23 /: f23[g23[x_]] := outer23[x]");
    assert_eval_eq("f23[g23[1]]", "outer23[1]", 0);
    ev("Block[{g23}, g23 /: f23[g23[x_]] := inner23[x]]");
    assert_eval_eq("f23[g23[1]]", "outer23[1]", 0);   /* restored */
}

/* ---- Clear / Remove drop upvalues ---- */
static void test_clear_removes_upvalues(void) {
    ev("g24 /: f24[g24[x_]] := h24[x]");
    assert_eval_eq("f24[g24[3]]", "h24[3]", 0);
    ev("Clear[g24]");
    assert_eval_eq("UpValues[g24]", "{}", 0);
    assert_eval_eq("f24[g24[3]]", "f24[g24[3]]", 0);
}

/* ---- efficiency gate counter stays consistent ---- */
static void test_upvalue_count_consistency(void) {
    extern size_t symtab_up_value_count;
    size_t before = symtab_up_value_count;
    ev("effc25 /: wrap25[effc25[x_]] := x");
    assert(symtab_up_value_count == before + 1);
    ev("Clear[effc25]");
    assert(symtab_up_value_count == before);
}

/* ---- regression: Condition on a delayed RHS is lifted onto the LHS ---- */
static void test_condition_lift_upsetdelayed(void) {
    ev("fz26[gz26[x_]] ^:= yesz26[x] /; x > 0");
    assert_eval_eq("fz26[gz26[4]]", "yesz26[4]", 0);        /* guard holds */
    assert_eval_eq("fz26[gz26[-3]]", "fz26[gz26[-3]]", 0);  /* guard fails -> unevaluated */
}

static void test_condition_lift_tagsetdelayed(void) {
    ev("gc27 /: fc27[gc27[x_]] := pos27[x] /; x > 0");
    assert_eval_eq("fc27[gc27[5]]", "pos27[5]", 0);
    assert_eval_eq("fc27[gc27[-1]]", "fc27[gc27[-1]]", 0);
}

/* ---- regression: list setters honour Protected ---- */
static void test_protected_list_setter(void) {
    ev("Quiet[DownValues[Sin] = {HoldPattern[Sin[7]] :> 100}]");
    assert_eval_eq("Sin[7]", "Sin[7]", 0);                 /* builtin untouched */
    ev("Quiet[OwnValues[Pi] = {HoldPattern[Pi] :> 3}]");
    assert_eval_eq("Pi", "Pi", 0);                         /* constant untouched */
}

/* ---- regression: UpSet on pattern args does not pollute Blank / pattern heads ---- */
static void test_upset_no_blank_pollution(void) {
    ev("gn28 /: ffn28[x_, gn28] := hitn28");
    assert_eval_eq("ffn28[anything, gn28]", "hitn28", 0);  /* installed on the real symbol */
    assert_eval_eq("Length[UpValues[Blank]]", "0", 0);     /* NOT on Blank */
    /* a pattern-test element must not try to install on PatternTest */
    ev("gpt28 /: fpt28[x_?NumberQ, gpt28] := okpt28");
    assert_eval_eq("fpt28[3, gpt28]", "okpt28", 0);
}

static void test_upset_nosym(void) {
    /* f[x_] ^= rhs has no operand symbol -> message + left unevaluated */
    ev("Quiet[gnos28[x_] ^= 1]");
    assert_eval_eq("gnos28[5]", "gnos28[5]", 0);
}

/* ---- regression: upvalues on an Orderless head see WMA's sorted order ---- */
static void test_orderless_before_upvalue(void) {
    ev("rga29 /: rga29[x_] + rga29[y_] := rpl29[x, y]");
    /* args sorted (aaa < zzz) before the upvalue binds, so x=aaa, y=zzz --
     * independent of the written order (both match WMA's canonical binding). */
    assert_eval_eq("rga29[zzz] + rga29[aaa]", "rpl29[aaa, zzz]", 0);
    assert_eval_eq("rga29[aaa] + rga29[zzz]", "rpl29[aaa, zzz]", 0);
    /* a literal-head upvalue must fire regardless of the operand's written order
     * (regression: the first-arg-head dispatch filter is unsound for Orderless) */
    ev("rx29 /: rx29 + y_ := rf29[y]");
    assert_eval_eq("rx29 + 5", "rf29[5]", 0);
    assert_eval_eq("5 + rx29", "rf29[5]", 0);
}

/* ---- regression: Unset must not CORRUPT a non-matching guarded rule ----
 * (pattern_alpha_normalize used to rename the stored pattern's variables in place
 * via a refcount-aliased copy; a non-matching `=.` left a half-renamed rule.) */
static void test_tagunset_condition_no_corruption(void) {
    ev("cg31 /: cf31[cg31[x_]] := rr31[x] /; x > 0");
    /* unconditioned `=.` does NOT match the conditioned rule (WMA) -> rule kept,
     * and crucially NOT corrupted: it must still fire with x bound. */
    ev("cg31 /: cf31[cg31[x_]] =.");
    assert_eval_eq("cf31[cg31[3]]", "rr31[3]", 0);
    /* removal WITH the matching condition does remove it */
    ev("cg31 /: cf31[cg31[x_]] /; x > 0 =.");
    assert_eval_eq("cf31[cg31[3]]", "cf31[cg31[3]]", 0);
}

static void test_downvalue_unset_condition_no_corruption(void) {
    ev("dv31[x_] := rr31b[x] /; x > 0");
    ev("dv31[x_] =.");               /* non-matching unset must not corrupt */
    assert_eval_eq("dv31[3]", "rr31b[3]", 0);
}

/* ---- regression: a downvalue on the operand masks the upvalue (WMA) ---- */
static void test_operand_downvalue_masks_upvalue(void) {
    ev("md30 /: pd30[md30[x_]] := up30[x]");
    ev("md30[x_] := down30[x]");
    /* md30[5] rewrites to down30[5] during arg eval, so pd30 never sees md30 */
    assert_eval_eq("pd30[md30[5]]", "pd30[down30[5]]", 0);
}

/* ---- Flat-head leftover matching: a Flat (Plus/Times/...) rule fires on a
 * longer operand sequence than the pattern specifies, rewriting a subset and
 * keeping the rest (Mathematica semantics; added when the review of PR #87
 * flagged that the headline mod example only worked at exact arity). ---- */
static void test_flat_upvalue_longer_sum(void) {
    ev("mod70 /: mod70[a_, p_] + mod70[b_, p_] := mod70[Mod[a + b, p], p]");
    assert_eval_eq("mod70[1,5] + mod70[2,5]", "mod70[3,5]", 0);               /* exact arity */
    assert_eval_eq("mod70[1,5] + mod70[2,5] + mod70[3,5]", "mod70[1,5]", 0);  /* leftover -> fixed point */
    assert_eval_eq("mod70[1,5] + mod70[2,5] + z70", "mod70[3,5] + z70", 0);   /* non-matching leftover kept */
    ev("d70 /: d70[a_] + d70[b_] := d70[a + b]");
    assert_eval_eq("d70[1] + d70[2] + d70[3]", "d70[6]", 0);
    assert_eval_eq("d70[1] + d70[2] + d70[3] + d70[4]", "d70[10]", 0);
}

static void test_flat_downvalue_leftover(void) {
    /* a user Flat+Orderless head: DownValue (not UpValue) fires with leftover */
    ev("SetAttributes[g71, {Flat, Orderless}]");
    ev("g71[k71[a_], k71[b_]] := k71[a + b]");
    assert_eval_eq("g71[k71[1], k71[2], m71]", "g71[k71[3], m71]", 0);
}

static void test_flat_ordered_prefix_only(void) {
    /* non-Orderless Flat: the trailing-only leftover catcher means matched
     * elements must form a PREFIX run; an interior match does not fire. */
    ev("SetAttributes[hh72, Flat]");
    ev("hh72[u72[a_], u72[b_]] := u72[a + b]");
    assert_eval_eq("hh72[u72[1], u72[2], w72]", "hh72[u72[3], w72]", 0);        /* prefix fires */
    assert_eval_eq("hh72[w72, u72[1], u72[2]]", "hh72[w72, u72[1], u72[2]]", 0);/* interior: unevaluated */
}

static void test_flat_group_absorption_not_augment(void) {
    /* MUST stay group-absorption (NOT leftover): a trailing blank swallows the
     * surplus as ONE group, so a_ binds x73+y73, not x73 with y73 left over. */
    ev("pp73 /: Plus[pp73, a_] := lab73[a]");
    assert_eval_eq("pp73 + x73", "lab73[x73]", 0);
    assert_eval_eq("pp73 + x73 + y73", "lab73[x73 + y73]", 0);
}

static void test_upset_locked(void) {
    /* #2 from the review: UpSet must honour Locked, not just Protected. */
    ev("SetAttributes[lk74, Locked]");
    ev("foo74[lk74] ^= 3");                 /* refused: UpSet::write (Locked) */
    assert_eval_eq("UpValues[lk74]", "{}", 0);
}

static void test_flat_matchq_whole_expression(void) {
    /* Leftover matching is a rule-application feature, NOT a matcher change:
     * MatchQ must stay a WHOLE-expression test (matches WL). */
    assert_eval_eq("MatchQ[modq[1,5] + modq[2,5] + xq, modq[a_,p_] + modq[b_,p_]]", "False", 0);
    assert_eval_eq("MatchQ[aq + bq + cq, x_ + y_]", "True", 0);   /* group absorption */
    assert_eval_eq("MatchQ[aq + bq, x_ + y_]", "True", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_upset_basic);
    TEST(test_upset_multi_symbol);
    TEST(test_upset_immediate_rhs);
    TEST(test_upsetdelayed);
    TEST(test_upsetdelayed_arg_head);
    TEST(test_tagsetdelayed_upvalue);
    TEST(test_tagset_downvalue_redundant);
    TEST(test_tagset_ownvalue_redundant);
    TEST(test_tagset_blank_head);
    TEST(test_tagset_not_found);
    TEST(test_modular_arithmetic);
    TEST(test_tagunset);
    TEST(test_plain_unset_keeps_upvalue);
    TEST(test_upvalue_before_downvalue);
    TEST(test_allheads_hold);
    TEST(test_upvalues_reader);
    TEST(test_upvalues_setter);
    TEST(test_upvalues_copy);
    TEST(test_upvalues_reorder);
    TEST(test_downvalues_setter);
    TEST(test_ownvalues_setter);
    TEST(test_definition_inert_fullform);
    TEST(test_block_restores_upvalues);
    TEST(test_clear_removes_upvalues);
    TEST(test_upvalue_count_consistency);

    /* regressions found during stress testing */
    TEST(test_condition_lift_upsetdelayed);
    TEST(test_condition_lift_tagsetdelayed);
    TEST(test_protected_list_setter);
    TEST(test_upset_no_blank_pollution);
    TEST(test_upset_nosym);
    TEST(test_orderless_before_upvalue);
    TEST(test_tagunset_condition_no_corruption);
    TEST(test_downvalue_unset_condition_no_corruption);
    TEST(test_operand_downvalue_masks_upvalue);

    /* Flat-head leftover matching + Locked (PR #87 review follow-ups) */
    TEST(test_flat_upvalue_longer_sum);
    TEST(test_flat_downvalue_leftover);
    TEST(test_flat_ordered_prefix_only);
    TEST(test_flat_group_absorption_not_augment);
    TEST(test_upset_locked);
    TEST(test_flat_matchq_whole_expression);

    printf("All UpValues tests passed.\n");
    return 0;
}
