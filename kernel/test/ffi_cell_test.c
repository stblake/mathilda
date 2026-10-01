#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "ffi/mathilda_ffi.h"

static char g_lines[256][4096];
static int  g_n = 0;
static void sink(void* ctx, const char* line) {
    (void)ctx;
    if (g_n < 256) { strncpy(g_lines[g_n], line, 4095); g_lines[g_n][4095]=0; g_n++; }
}
static void reset(void){ g_n = 0; }

/* crude: extract "type":"X" from a line */
static const char* linetype(const char* l){
    static char t[32]; const char* p = strstr(l, "\"type\":\""); if(!p) return "?";
    p += 8; int i=0; while(*p && *p!='"' && i<31) t[i++]=*p++; t[i]=0; return t;
}
static int fails = 0;
static void expect_kinds(const char* name, const char* expect){
    char got[512]=""; for(int i=0;i<g_n;i++){ if(i)strcat(got,","); strcat(got, linetype(g_lines[i])); }
    int ok = strcmp(got, expect)==0;
    printf("%s  %s\n", ok?"PASS":"FAIL", name);
    if(!ok){ printf("      expect: %s\n      got:    %s\n", expect, got); fails++; }
}
static int any_contains(const char* sub){ for(int i=0;i<g_n;i++) if(strstr(g_lines[i],sub)) return 1; return 0; }
static void check(const char* name, int cond){ printf("%s  %s\n", cond?"PASS":"FAIL", name); if(!cond) fails++; }

int main(int argc, char** argv){
    /* The dir containing internal/ (init.m): argv[1], else $MATHILDA_HOME,
     * else "src" (run from the repo root). */
    const char* home = (argc > 1) ? argv[1]
                     : (getenv("MATHILDA_HOME") ? getenv("MATHILDA_HOME") : "src");
    mathilda_ffi_set_home(home);
    mathilda_ffi_init();
    printf("version: %s\n", mathilda_ffi_version());

    /* ---- is_complete ---- */
    check("is_complete: balanced", mathilda_ffi_is_complete("f[x] + {1,2}")==1);
    check("is_complete: open bracket", mathilda_ffi_is_complete("f[1,")==0);
    check("is_complete: open string", mathilda_ffi_is_complete("\"abc")==0);
    check("is_complete: open comment", mathilda_ffi_is_complete("(* hi")==0);
    check("is_complete: bracket in string ok", mathilda_ffi_is_complete("\"f[\"")==1);

    /* ---- complete ---- */
    char* c = mathilda_ffi_complete("Sin");
    printf("complete(\"Sin\") = %s\n", c);
    check("complete: has Sin", strstr(c,"\"Sin\"")!=NULL);
    check("complete: has Sinh", strstr(c,"\"Sinh\"")!=NULL);
    check("complete: has SinIntegral", strstr(c,"\"SinIntegral\"")!=NULL);
    check("complete: sorted (Sin before Sinh)", strstr(c,"\"Sin\"")<strstr(c,"\"Sinh\""));
    mathilda_ffi_free(c);

    /* ---- eval_cell: mirror check_pipe_protocol cases ---- */
    reset(); mathilda_ffi_eval_cell("a = 1\nb = 2\na + b", 1, 1, sink, NULL);
    expect_kinds("cell: one stmt per line, every result shown",
                 "line,expr,line,expr,line,expr,done");

    reset(); mathilda_ffi_eval_cell("Print[\"hello\"]; Print[\"world\"]; 1/0", 2, 1, sink, NULL);
    expect_kinds("cell: Print->stream before result; 1/0 message",
                 "line,stream,line,stream,line,message,expr,done");
    { char j[256]=""; for(int i=0;i<g_n;i++){ const char*t=g_lines[i]; if(strstr(t,"\"type\":\"stream\"")){ const char*q=strstr(t,"\"text\":\""); if(q){q+=8; while(*q&&*q!='"'){ if(q[0]=='\\'&&q[1]=='n'){strcat(j,"\n");q+=2;} else {size_t L=strlen(j); j[L]=*q++; j[L+1]=0;} } } } }
      check("cell: joined stream text == hello\\nworld\\n", strcmp(j,"hello\nworld\n")==0); }
    check("cell: Power::infy message", any_contains("Power::infy"));

    reset(); mathilda_ffi_eval_cell("x = 5;", 3, 1, sink, NULL);
    expect_kinds("cell: x=5; sends only its line", "line,done");

    reset(); mathilda_ffi_eval_cell("x = 5; y = 6", 4, 1, sink, NULL);
    expect_kinds("cell: x=5; y=6 shows only y", "line,line,expr,done");

    reset(); mathilda_ffi_eval_cell("f[x_] :=\n  x^2\nf[3]", 5, 1, sink, NULL);
    expect_kinds("cell: stmt continues over line break", "line,line,expr,done");
    check("cell: f[3] -> 9", any_contains("\"payload\":\"9\""));

    reset(); mathilda_ffi_eval_cell("1 + 1\nf[1,", 6, 1, sink, NULL);
    expect_kinds("cell: syntax error anywhere evaluates nothing", "error,done");

    reset(); mathilda_ffi_eval_cell("(* just a comment *)", 7, 1, sink, NULL);
    expect_kinds("cell: comment-only cell is just done", "done");

    reset(); mathilda_ffi_eval_cell("Print[1]\nSin[1, 2]\n3", 8, 1, sink, NULL);
    expect_kinds("cell: per-statement order", "line,stream,line,message,expr,line,expr,done");

    reset(); mathilda_ffi_eval_cell("?Sin", 9, 1, sink, NULL);
    expect_kinds("cell: ?Sin is a usage message", "line,usage,done");

    /* ---- history: %, %%, Out[n] ---- */
    reset(); mathilda_ffi_eval_cell("11 + 11", 1, 1, sink, NULL); expect_kinds("hist seed", "line,expr,done");
    reset(); mathilda_ffi_eval_cell("% + 1", 2, 1, sink, NULL);
    check("cell: % is previous result (23)", any_contains("\"payload\":\"23\""));
    reset(); mathilda_ffi_eval_cell("%%", 3, 1, sink, NULL);
    check("cell: %% two back (22)", any_contains("\"payload\":\"22\""));
    /* Out[seed_line]: find the $Line the seed took (session counter is cumulative). */
    reset(); mathilda_ffi_eval_cell("11 + 11", 5, 1, sink, NULL);
    int seedline=-1; for(int i=0;i<g_n;i++){ const char*q=strstr(g_lines[i],"\"line\":"); if(q&&strstr(g_lines[i],"\"type\":\"line\"")){ seedline=atoi(q+7); } }
    { char buf[64]; snprintf(buf,sizeof buf,"Out[%d]",seedline); reset(); mathilda_ffi_eval_cell(buf,6,1,sink,NULL);
      check("cell: Out[seedline] == 22", any_contains("\"payload\":\"22\"")); }

    /* ---- plain mode (cell=0): Null shown, no stream ---- */
    reset(); mathilda_ffi_eval_cell("Print[\"x\"]; 7", 10, 0, sink, NULL);
    expect_kinds("plain: single expr, no stream/line", "expr,done");

    /* ---- latex present on a fraction ---- */
    reset(); mathilda_ffi_eval_cell("1/2", 11, 1, sink, NULL);
    check("cell: latex field present for 1/2", any_contains("\"latex\":"));

    printf("\n%s (%d failure%s)\n", fails?"FAILURES":"ALL PASS", fails, fails==1?"":"s");
    return fails?1:0;
}
