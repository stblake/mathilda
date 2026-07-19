/* Tests for ListConvolve — discrete convolution.
 *
 * Golden vectors from the Wolfram Language reference cover: symbolic 1-D
 * convolution with every overhang form (default, single index, {kL,kR}),
 * constant / cyclic-list / empty padding, generalized g/h in place of
 * Times/Plus (including empty-pad term dropping and the h argument order),
 * multidimensional kernels/data with per-axis overhangs, exact-rational data,
 * the FFT numeric fast path (machine + MPFR, 1-D and 2-D) checked against a
 * definitional reference, argument diagnostics, and attributes. */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* Evaluate `input` and `expected`; assert structural equality after evaluation
 * (both are canonicalised by the evaluator, so operand order is normalised). */
static void chk_eq(const char* input, const char* expected) {
    Expr* ie = parse_expression(input);
    ASSERT_MSG(ie != NULL, "parse failed: %s", input);
    Expr* iv = evaluate(ie);
    expr_free(ie);
    Expr* ee = parse_expression(expected);
    ASSERT_MSG(ee != NULL, "parse failed: %s", expected);
    Expr* ev = evaluate(ee);
    expr_free(ee);
    if (!expr_eq(iv, ev)) {
        char* a = expr_to_string(iv);
        char* b = expr_to_string(ev);
        ASSERT_MSG(0, "%s\n  expected: %s\n  actual:   %s", input, b, a);
        free(a); free(b);
    }
    expr_free(iv);
    expr_free(ev);
}

/* Assert a numeric expression evaluates to a real <= tol in magnitude. */
static void chk_small(const char* input, double tol) {
    Expr* e = parse_expression(input);
    ASSERT_MSG(e != NULL, "parse failed: %s", input);
    Expr* r = evaluate(e);
    expr_free(e);
    double v = NAN;
    if (r->type == EXPR_REAL) v = r->data.real;
    else if (r->type == EXPR_INTEGER) v = (double)r->data.integer;
#ifdef USE_MPFR
    else if (r->type == EXPR_MPFR) v = mpfr_get_d(r->data.mpfr, MPFR_RNDN);
#endif
    ASSERT_MSG(!isnan(v) && fabs(v) <= tol, "%s: expected <= %g, got %g",
               input, tol, v);
    expr_free(r);
}

/* Assert `input` is returned unevaluated with its head intact — i.e. it prints
 * as `head[...]`. Unlike chk_eq(X, X) (a tautology, since both sides evaluate
 * the same way), this genuinely distinguishes "left unevaluated" from "reduced
 * to a result". */
static void chk_stays(const char* input, const char* head) {
    Expr* e = parse_expression(input);
    ASSERT_MSG(e != NULL, "parse failed: %s", input);
    Expr* r = evaluate(e);
    expr_free(e);
    char* s = expr_to_string(r);
    size_t hl = strlen(head);
    ASSERT_MSG(strncmp(s, head, hl) == 0 && s[hl] == '[',
               "%s: expected unevaluated %s[...], got %s", input, head, s);
    free(s);
    expr_free(r);
}

/* ---- 1-D symbolic, overhangs ---------------------------------------- */

static void test_basic(void) {
    chk_eq("ListConvolve[{x,y},{a,b,c,d,e,f}]",
           "{b x+a y,c x+b y,d x+c y,e x+d y,f x+e y}");
    chk_eq("ListConvolve[{x,y,z},{1,2,3,4,5,6}]",
           "{3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,6 x+5 y+4 z}");
    /* Default output is not the same length as the input. */
    chk_eq("ListConvolve[{1,1},{1,2,3,4,5}]", "{3,5,7,9}");
    chk_eq("ListConvolve[{1,1}/2,{a,b,c,d,e}]",
           "{a/2+b/2,b/2+c/2,c/2+d/2,d/2+e/2}");
}

static void test_overhangs(void) {
    chk_eq("ListConvolve[{x,y},{a,b,c,d,e,f},1]",
           "{a x+f y,b x+a y,c x+b y,d x+c y,e x+d y,f x+e y}");
    chk_eq("ListConvolve[{x,y},{a,b,c,d,e,f},2]",
           "{b x+a y,c x+b y,d x+c y,e x+d y,f x+e y,a x+f y}");
    chk_eq("ListConvolve[{x,y,z},{1,2,3,4,5,6},1]",
           "{x+6 y+5 z,2 x+y+6 z,3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,6 x+5 y+4 z}");
    chk_eq("ListConvolve[{x,y,z},{1,2,3,4,5,6},-1]",
           "{3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,6 x+5 y+4 z,x+6 y+5 z,2 x+y+6 z}");
    chk_eq("ListConvolve[{x,y,z},{1,2,3,4,5,6},{1,-1}]",
           "{x+6 y+5 z,2 x+y+6 z,3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,"
           "6 x+5 y+4 z,x+6 y+5 z,2 x+y+6 z}");
    /* Make the output the same length as the input. */
    chk_eq("ListConvolve[{1,1},{1,2,3,4,5},1]", "{6,3,5,7,9}");
}

/* ---- padding -------------------------------------------------------- */

static void test_padding(void) {
    chk_eq("ListConvolve[{x,y},{a,b,c,d,e,f},1,zzz]",
           "{a x+y zzz,b x+a y,c x+b y,d x+c y,e x+d y,f x+e y}");
    chk_eq("ListConvolve[{x,y,z},{1,2,3,4,5},1,aa]",
           "{x+aa y+aa z,2 x+y+aa z,3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z}");
    chk_eq("ListConvolve[{x,y,z},{1,2,3,4,5},1,{aa,bb}]",
           "{x+bb y+aa z,2 x+y+bb z,3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z}");
    chk_eq("ListConvolve[{x,y,z},{1,2,3,4,5},1,{aa,bb,cc}]",
           "{x+cc y+bb z,2 x+y+cc z,3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z}");
    chk_eq("ListConvolve[{x,y,z},{1,2,3},{1,-1},{aa,bb,cc}]",
           "{x+cc y+bb z,2 x+y+cc z,3 x+2 y+z,aa x+3 y+2 z,bb x+aa y+3 z}");
}

/* ---- generalized g / h --------------------------------------------- */

static void test_generalized(void) {
    chk_eq("ListConvolve[{x,y},{1,2,3,4},1,0,f]",
           "{f[x,1]+f[y,0],f[x,2]+f[y,1],f[x,3]+f[y,2],f[x,4]+f[y,3]}");
    chk_eq("ListConvolve[{x,y},{1,2,3,4},1,0,f,g]",
           "{g[f[y,0],f[x,1]],g[f[y,1],f[x,2]],g[f[y,2],f[x,3]],g[f[y,3],f[x,4]]}");
    chk_eq("ListConvolve[{x,y},{1,2,3,4},1,0,f,List]",
           "{{f[y,0],f[x,1]},{f[y,1],f[x,2]},{f[y,2],f[x,3]},{f[y,3],f[x,4]}}");
    /* Empty padding drops the missing list factor (single-argument g term). */
    chk_eq("ListConvolve[{x,y,z},{1,2,3},1,{},f,g]",
           "{g[f[z],f[y],f[x,1]],g[f[z],f[y,1],f[x,2]],g[f[z,1],f[y,2],f[x,3]]}");
    chk_eq("ListConvolve[{x,y,z},{1,2,3},-1,{},f,g]",
           "{g[f[z,1],f[y,2],f[x,3]],g[f[z,2],f[y,3],f[x]],g[f[z,3],f[y],f[x]]}");
    chk_eq("ListConvolve[{x,y,z},{1,2,3},{1,-1},{},f,g]",
           "{g[f[z],f[y],f[x,1]],g[f[z],f[y,1],f[x,2]],g[f[z,1],f[y,2],f[x,3]],"
           "g[f[z,2],f[y,3],f[x]],g[f[z,3],f[y],f[x]]}");
}

/* ---- multidimensional ---------------------------------------------- */

static void test_multidim(void) {
    chk_eq("ListConvolve[{{1,1},{1,1}},{{a,b,c},{d,e,f},{g,h,i}}]",
           "{{a+b+d+e,b+c+e+f},{d+e+g+h,e+f+h+i}}");
    chk_eq("ListConvolve[{{1,1},{1,1}},{{1,2,3},{4,5,6},{7,8,9}},1]",
           "{{20,18,22},{14,12,16},{26,24,28}}");
    chk_eq("ListConvolve[{{a,b},{c,d}},{{1,2,3},{4,5,6},{7,8,9}},1]",
           "{{a+3 b+7 c+9 d,2 a+b+8 c+7 d,3 a+2 b+9 c+8 d},"
           "{4 a+6 b+c+3 d,5 a+4 b+2 c+d,6 a+5 b+3 c+2 d},"
           "{7 a+9 b+4 c+6 d,8 a+7 b+5 c+4 d,9 a+8 b+6 c+5 d}}");
    /* Per-axis overhangs {{kL1,kL2},{kR1,kR2}}. */
    chk_eq("ListConvolve[{{a,b},{c,d}},{{1,2,3},{4,5,6},{7,8,9}},{{1,1},{2,2}}]",
           "{{a+3 b+7 c+9 d,2 a+b+8 c+7 d,3 a+2 b+9 c+8 d,a+3 b+7 c+9 d},"
           "{4 a+6 b+c+3 d,5 a+4 b+2 c+d,6 a+5 b+3 c+2 d,4 a+6 b+c+3 d},"
           "{7 a+9 b+4 c+6 d,8 a+7 b+5 c+4 d,9 a+8 b+6 c+5 d,7 a+9 b+4 c+6 d},"
           "{a+3 b+7 c+9 d,2 a+b+8 c+7 d,3 a+2 b+9 c+8 d,a+3 b+7 c+9 d}}");
}

/* ---- exact numeric (direct path stays exact) ----------------------- */

static void test_exact(void) {
    chk_eq("ListConvolve[Reverse[Differences[Range[10]^2]], 1/Range[20]]",
           "ListCorrelate[Differences[Range[10]^2], 1/Range[20]]");
    chk_eq("ListConvolve[{1,2,3},{1/2,1/3,1/4,1/5,1/6}]",
           "{ListConvolve[{1,2,3},{1/2,1/3,1/4,1/5,1/6}][[1]],"
           " ListConvolve[{1,2,3},{1/2,1/3,1/4,1/5,1/6}][[2]],"
           " ListConvolve[{1,2,3},{1/2,1/3,1/4,1/5,1/6}][[3]]}");
    /* Concrete exact values: Sum_r ker[r] list[s-r]. */
    chk_eq("ListConvolve[{1,2,3},{10,20,30,40,50}]",
           "{3*10+2*20+1*30,3*20+2*30+1*40,3*30+2*40+1*50}");
}

/* ---- FFT numeric fast path vs definitional reference ---------------- */

static void test_fft_machine(void) {
    /* work = L*m = 171*30 >> threshold => machine FFT path. */
    chk_small(
        "SeedRandom[11]; ker=RandomReal[{-1,1},30]; lst=RandomReal[{-1,1},200];"
        "m=Length[ker]; n=Length[lst]; L=n-m+1;"
        "ref=Table[Total[Reverse[ker]*lst[[t;;t+m-1]]],{t,1,L}];"
        "Max[Abs[ListConvolve[ker,lst]-ref]]",
        1e-9);
    /* 2-D machine FFT (8x8 kernel, 40x40 data). */
    chk_small(
        "SeedRandom[12]; ker=RandomReal[{-1,1},{8,8}]; dat=RandomReal[{-1,1},{40,40}];"
        "m1=8;m2=8;L1=33;L2=33; rk=Reverse[Reverse[ker],2];"
        "ref=Table[Total[Flatten[rk*dat[[t1;;t1+m1-1,t2;;t2+m2-1]]]],{t1,1,L1},{t2,1,L2}];"
        "Max[Abs[Flatten[ListConvolve[ker,dat]-ref]]]",
        1e-9);
}

#ifdef USE_MPFR
static void test_fft_mpfr(void) {
    /* Genuine 30-digit inputs (exact rationals -> N[...,30]) exercise the MPFR
     * FFT and must retain high precision. */
    chk_small(
        "ker=N[Table[1/(r^2+1),{r,1,80}],30]; lst=N[Table[1/(s+3),{s,1,300}],30];"
        "m=80;L=221; ref=Table[Total[Reverse[ker]*lst[[t;;t+m-1]]],{t,1,L}];"
        "Max[Abs[ListConvolve[ker,lst]-ref]]",
        1e-25);
    /* Small MPFR case (below the FFT threshold) exercises the direct numeric
     * loop and must likewise retain 30-digit precision. */
    chk_small(
        "ker=N[Table[1/(r+1),{r,1,3}],30]; lst=N[Table[1/(s^2+1),{s,1,8}],30];"
        "m=3;L=6; ref=Table[Total[Reverse[ker]*lst[[t;;t+m-1]]],{t,1,L}];"
        "Max[Abs[ListConvolve[ker,lst]-ref]]",
        1e-27);
}
#endif

/* ---- kernel/list size relationships -------------------------------- */

static void test_kernel_sizes(void) {
    /* Single-element kernel scales the list (m = 1). */
    chk_eq("ListConvolve[{k},{a,b,c}]", "{a k,b k,c k}");
    chk_eq("ListConvolve[{3},{a,b,c,d}]", "{3 a,3 b,3 c,3 d}");
    /* Kernel the same length as the list: default overhang gives one output,
     * with the h-terms ordered by ascending list index. */
    chk_eq("ListConvolve[{x,y,z},{a,b,c}]", "{c x+b y+a z}");
    chk_eq("ListConvolve[{x,y},{a,b}]", "{b x+a y}");
}

/* ListConvolve[k,l] == ListCorrelate[Reverse[k],l] (default overhang); with an
 * explicit overhang the {kL,kR} convention negates in place under reversal. */
static void test_reverse_identity(void) {
    chk_eq("ListConvolve[{x,y,z},{a,b,c,d,e}]",
           "ListCorrelate[Reverse[{x,y,z}],{a,b,c,d,e}]");
    chk_eq("ListConvolve[{p,q,r,s},{1,2,3,4,5,6,7},{1,-1}]",
           "ListCorrelate[Reverse[{p,q,r,s}],{1,2,3,4,5,6,7},{-1,1}]");
}

/* ---- optional level argument (7th) --------------------------------- */

static void test_lev_argument(void) {
    /* lev == rank(ker) is the only supported value and is a no-op. */
    chk_eq("ListConvolve[{x,y},{a,b,c,d},1,p,Times,Plus,1]",
           "ListConvolve[{x,y},{a,b,c,d},1,p]");
    chk_eq("ListConvolve[{{a,b},{c,d}},{{1,2,3},{4,5,6},{7,8,9}},1,0,Times,Plus,2]",
           "ListConvolve[{{a,b},{c,d}},{{1,2,3},{4,5,6},{7,8,9}},1,0]");
    /* lev != rank(ker) is unsupported: the call stays unevaluated. */
    chk_stays("ListConvolve[{x,y},{a,b,c,d},1,p,Times,Plus,2]", "ListConvolve");
}

/* ---- numeric direct path (below the FFT threshold) ----------------- */

static void test_machine_direct(void) {
    /* Small Real input takes the tight O(L*m) numeric loop, not the FFT;
     * checked against the definitional reference Sum_r Reverse[ker]_r lst. */
    chk_small(
        "ker={0.5,-1.5,2.0}; lst={1.0,2.0,3.0,4.0,5.0,6.0};"
        "m=3; L=4; ref=Table[Total[Reverse[ker]*lst[[t;;t+m-1]]],{t,1,L}];"
        "Max[Abs[ListConvolve[ker,lst]-ref]]",
        1e-12);
    /* Complex machine data through the same loop. */
    chk_small(
        "ker={1.0+1.0 I, 2.0-1.0 I}; lst={0.5, 1.0+2.0 I, -1.0, 3.0 I, 2.0};"
        "m=2; L=4; ref=Table[Total[Reverse[ker]*lst[[t;;t+m-1]]],{t,1,L}];"
        "Max[Abs[ListConvolve[ker,lst]-ref]]",
        1e-12);
}

/* ---- rank-3 (exact engine, all-ones kernel) ------------------------ */

static void test_rank3(void) {
    /* 2x2x2 all-ones kernel over a 3x3x3 array: each output is the sum over a
     * 2x2x2 window (reversal is irrelevant for an all-ones kernel). Exact. */
    chk_small(
        "dat=Table[100 i+10 j+k,{i,1,3},{j,1,3},{k,1,3}];"
        "ker=ConstantArray[1,{2,2,2}];"
        "ref=Table[Sum[dat[[t1+i,t2+j,t3+k]],{i,0,1},{j,0,1},{k,0,1}],"
        "{t1,1,2},{t2,1,2},{t3,1,2}];"
        "Max[Abs[Flatten[ListConvolve[ker,dat]-ref]]]",
        0.0);
}

/* ---- diagnostics + attributes -------------------------------------- */

static void test_errors_attrs(void) {
    /* Wrong arity: stays unevaluated. */
    chk_stays("ListConvolve[{1,2}]", "ListConvolve");
    chk_stays("ListConvolve[a,b,c,d,e,f,g,h]", "ListConvolve");  /* 8 args */
    /* Non-list arguments: unevaluated. */
    chk_stays("ListConvolve[x,y]", "ListConvolve");
    /* Kernel and list of different ranks. */
    chk_stays("ListConvolve[{x,y},{{a,b},{c,d}}]", "ListConvolve");
    /* Zero overhang index is invalid. */
    chk_stays("ListConvolve[{x,y},{a,b,c},0]", "ListConvolve");
    /* An overhang that empties the output window (L = n - KL + KR = 0). */
    chk_stays("ListConvolve[{v,w,x,y,z},{a,b},{3,1}]", "ListConvolve");
    /* Empty kernel. */
    chk_stays("ListConvolve[{},{a,b,c}]", "ListConvolve");
    /* Pad list of the wrong rank. */
    chk_stays("ListConvolve[{{a,b},{c,d}},{{1,2},{3,4}},1,{9,9}]", "ListConvolve");
    ASSERT(symtab_get_def("ListConvolve")->attributes & ATTR_PROTECTED);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_basic);
    TEST(test_overhangs);
    TEST(test_padding);
    TEST(test_generalized);
    TEST(test_multidim);
    TEST(test_exact);
    TEST(test_kernel_sizes);
    TEST(test_reverse_identity);
    TEST(test_lev_argument);
    TEST(test_machine_direct);
    TEST(test_rank3);
    TEST(test_fft_machine);
#ifdef USE_MPFR
    TEST(test_fft_mpfr);
#endif
    TEST(test_errors_attrs);

    printf("All ListConvolve tests passed.\n");
    return 0;
}
