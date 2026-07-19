/* Tests for ListCorrelate — discrete correlation.
 *
 * Golden vectors from the Wolfram Language reference cover: symbolic 1-D
 * correlation with every overhang form, the ListCorrelate/ListConvolve[Reverse]
 * identity, constant / empty padding, generalized g/h, multidimensional data,
 * exact-rational data, the FFT numeric fast path (machine + MPFR) checked
 * against a definitional reference, complex data, diagnostics, and attributes. */

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

/* Assert `input` is returned unevaluated with its head intact (prints as
 * `head[...]`). Unlike chk_eq(X, X) — a tautology — this genuinely checks the
 * call was left unevaluated rather than reduced. */
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
    chk_eq("ListCorrelate[{x,y},{a,b,c,d,e,f}]",
           "{a x+b y,b x+c y,c x+d y,d x+e y,e x+f y}");
    chk_eq("ListCorrelate[{x,y},{a,b,c,d,e,f},1]",
           "{a x+b y,b x+c y,c x+d y,d x+e y,e x+f y,f x+a y}");
    chk_eq("ListCorrelate[{x,y},{a,b,c,d,e,f},2]",
           "{f x+a y,a x+b y,b x+c y,c x+d y,d x+e y,e x+f y}");
}

/* ListCorrelate[ker,list] == ListConvolve[Reverse[ker],list] (default overhang;
 * the {kL,kR} conventions are negated between the two, so with explicit
 * overhangs the reversed kernel must use the negated {-kR,-kL}). */
static void test_reverse_identity(void) {
    chk_eq("ListCorrelate[{x,y,z},{a,b,c,d,e}]",
           "ListConvolve[Reverse[{x,y,z}],{a,b,c,d,e}]");
    chk_eq("ListCorrelate[{p,q,r,s},{1,2,3,4,5,6,7},{1,1}]",
           "ListConvolve[Reverse[{p,q,r,s}],{1,2,3,4,5,6,7},{-1,-1}]");
}

/* ---- padding + generalized g/h ------------------------------------- */

static void test_padding_generalized(void) {
    chk_eq("ListCorrelate[{x,y},{a,b,c,d,e,f},1,zzz]",
           "{a x+b y,b x+c y,c x+d y,d x+e y,e x+f y,f x+y zzz}");
    chk_eq("ListCorrelate[{x,y,z},{1,2,3,4,5},{1,-1},{x,y,z},f,g]",
           "{g[f[x,1],f[y,2],f[z,3]],g[f[x,2],f[y,3],f[z,4]],g[f[x,3],f[y,4],f[z,5]]}");
    chk_eq("ListCorrelate[{x,y,z},{1,2,3,4,5},{-1,1},0,f,g]",
           "{g[f[x,0],f[y,0],f[z,1]],g[f[x,0],f[y,1],f[z,2]],g[f[x,1],f[y,2],f[z,3]],"
           "g[f[x,2],f[y,3],f[z,4]],g[f[x,3],f[y,4],f[z,5]],g[f[x,4],f[y,5],f[z,0]],"
           "g[f[x,5],f[y,0],f[z,0]]}");
    chk_eq("ListCorrelate[{x,y,z},{1,2,3,4,5},{-1,1},{},f,g]",
           "{g[f[x],f[y],f[z,1]],g[f[x],f[y,1],f[z,2]],g[f[x,1],f[y,2],f[z,3]],"
           "g[f[x,2],f[y,3],f[z,4]],g[f[x,3],f[y,4],f[z,5]],g[f[x,4],f[y,5],f[z]],"
           "g[f[x,5],f[y],f[z]]}");
}

/* ---- multidimensional ---------------------------------------------- */

static void test_multidim(void) {
    chk_eq("ListCorrelate[{{1,1},{1,1}},{{a,b,c},{d,e,f},{g,h,i}}]",
           "{{a+b+d+e,b+c+e+f},{d+e+g+h,e+f+h+i}}");
    /* 2-D correlation with the {{1,0,1},{0,-4,0},{1,0,1}} Laplacian stencil,
     * no overhang, over a 4x4 array of distinct symbols a[i,j]. */
    chk_eq("ListCorrelate[{{1,0,1},{0,-4,0},{1,0,1}}, Array[a, {4,4}]]",
           "{{a[1,1]+a[1,3]+a[3,1]+a[3,3]-4 a[2,2],"
           "  a[1,2]+a[1,4]+a[3,2]+a[3,4]-4 a[2,3]},"
           " {a[2,1]+a[2,3]+a[4,1]+a[4,3]-4 a[3,2],"
           "  a[2,2]+a[2,4]+a[4,2]+a[4,4]-4 a[3,3]}}");
}

/* ---- exact-rational golden vector (request example) ----------------- */

static void test_exact(void) {
    chk_eq("ListCorrelate[Differences[Range[10]^2], 1/Range[20]]",
           "{52489/2520,40499/2520,124189/9240,64591/5544,531397/51480,371809/40040,"
           "55361/6552,86017/11088,5860303/816816,962251/144144,76494941/12252240,"
           "5910403/1007760}");
}

/* ---- FFT numeric fast path ----------------------------------------- */

static void test_fft(void) {
    chk_small(
        "SeedRandom[21]; ker=RandomReal[{-1,1},30]; lst=RandomReal[{-1,1},200];"
        "m=Length[ker]; n=Length[lst]; L=n-m+1;"
        "ref=Table[Total[ker*lst[[t;;t+m-1]]],{t,1,L}];"
        "Max[Abs[ListCorrelate[ker,lst]-ref]]",
        1e-9);
#ifdef USE_MPFR
    chk_small(
        "ker=N[Table[1/(r+2),{r,1,80}],30]; lst=N[Table[1/(s^2+1),{s,1,300}],30];"
        "m=80;L=221; ref=Table[Total[ker*lst[[t;;t+m-1]]],{t,1,L}];"
        "Max[Abs[ListCorrelate[ker,lst]-ref]]",
        1e-25);
#endif
}

/* Complex machine data. */
static void test_complex(void) {
    chk_small(
        "Max[Abs[ListCorrelate[{1,I},{0.5+0.2 I,0.1-0.3 I,2.0,1.0+1.0 I}]"
        " - {0.8+0.3 I,0.1+1.7 I,1.0+1.0 I}]]",
        1e-9);
}

/* ---- kernel/list size relationships -------------------------------- */

static void test_kernel_sizes(void) {
    /* Single-element kernel scales the list (m = 1). */
    chk_eq("ListCorrelate[{k},{a,b,c}]", "{a k,b k,c k}");
    chk_eq("ListCorrelate[{3},{a,b,c,d}]", "{3 a,3 b,3 c,3 d}");
    /* Kernel the same length as the list: default overhang gives one output,
     * with the h-terms ordered by ascending list index. */
    chk_eq("ListCorrelate[{x,y,z},{a,b,c}]", "{a x+b y+c z}");
    chk_eq("ListCorrelate[{x,y},{a,b}]", "{a x+b y}");
    /* Correlate/convolve reverse identity with an explicit overhang. */
    chk_eq("ListCorrelate[{x,y,z},{a,b,c,d,e},1]",
           "ListConvolve[Reverse[{x,y,z}],{a,b,c,d,e},-1]");
}

/* ---- optional level argument (7th) --------------------------------- */

static void test_lev_argument(void) {
    /* lev == rank(ker) is the only supported value and is a no-op. */
    chk_eq("ListCorrelate[{x,y},{a,b,c,d},1,p,Times,Plus,1]",
           "ListCorrelate[{x,y},{a,b,c,d},1,p]");
    chk_eq("ListCorrelate[{{a,b},{c,d}},{{1,2,3},{4,5,6},{7,8,9}},1,0,Times,Plus,2]",
           "ListCorrelate[{{a,b},{c,d}},{{1,2,3},{4,5,6},{7,8,9}},1,0]");
    /* lev != rank(ker) is unsupported: the call stays unevaluated. */
    chk_stays("ListCorrelate[{x,y},{a,b,c,d},1,p,Times,Plus,2]", "ListCorrelate");
}

/* ---- numeric direct path (below the FFT threshold) ----------------- */

static void test_machine_direct(void) {
    /* Small Real input takes the tight O(L*m) numeric loop, not the FFT;
     * checked against the definitional reference Sum_r ker_r lst. */
    chk_small(
        "ker={0.5,-1.5,2.0}; lst={1.0,2.0,3.0,4.0,5.0,6.0};"
        "m=3; L=4; ref=Table[Total[ker*lst[[t;;t+m-1]]],{t,1,L}];"
        "Max[Abs[ListCorrelate[ker,lst]-ref]]",
        1e-12);
}

/* ---- rank-3 (exact engine, all-ones kernel) ------------------------ */

static void test_rank3(void) {
    /* 2x2x2 all-ones kernel over a 3x3x3 array: each output is the sum over a
     * 2x2x2 window. Exact integer arithmetic (direct engine). */
    chk_small(
        "dat=Table[100 i+10 j+k,{i,1,3},{j,1,3},{k,1,3}];"
        "ker=ConstantArray[1,{2,2,2}];"
        "ref=Table[Sum[dat[[t1+i,t2+j,t3+k]],{i,0,1},{j,0,1},{k,0,1}],"
        "{t1,1,2},{t2,1,2},{t3,1,2}];"
        "Max[Abs[Flatten[ListCorrelate[ker,dat]-ref]]]",
        0.0);
}

/* ---- diagnostics + attributes -------------------------------------- */

static void test_errors_attrs(void) {
    /* Wrong arity: stays unevaluated. */
    chk_stays("ListCorrelate[{1,2}]", "ListCorrelate");
    chk_stays("ListCorrelate[a,b,c,d,e,f,g,h]", "ListCorrelate");  /* 8 args */
    /* Non-list arguments: unevaluated. */
    chk_stays("ListCorrelate[x,y]", "ListCorrelate");
    /* Kernel and list of different ranks. */
    chk_stays("ListCorrelate[{x,y},{{a,b},{c,d}}]", "ListCorrelate");
    /* Zero overhang index is invalid. */
    chk_stays("ListCorrelate[{x,y},{a,b,c},0]", "ListCorrelate");
    /* An overhang that empties the output window (L = n + KL - KR = 0). */
    chk_stays("ListCorrelate[{v,w,x,y,z},{a,b},{1,3}]", "ListCorrelate");
    /* Empty kernel. */
    chk_stays("ListCorrelate[{},{a,b,c}]", "ListCorrelate");
    ASSERT(symtab_get_def("ListCorrelate")->attributes & ATTR_PROTECTED);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_basic);
    TEST(test_reverse_identity);
    TEST(test_padding_generalized);
    TEST(test_multidim);
    TEST(test_exact);
    TEST(test_kernel_sizes);
    TEST(test_lev_argument);
    TEST(test_machine_direct);
    TEST(test_rank3);
    TEST(test_fft);
    TEST(test_complex);
    TEST(test_errors_attrs);

    printf("All ListCorrelate tests passed.\n");
    return 0;
}
