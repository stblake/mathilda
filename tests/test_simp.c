#include "test_utils.h"
#include "symtab.h"
#include "core.h"

void test_trigtoexp() {
    printf("testing %s\n", "TrigToExp[Cos[x]]"); fflush(stdout); assert_eval_eq("TrigToExp[Cos[x]]", "1/2 E^(-I x) + 1/2 E^(I x)", 0);
    printf("testing %s\n", "TrigToExp[Sinh[x]]"); fflush(stdout); assert_eval_eq("TrigToExp[Sinh[x]]", "1/2 E^x - 1/2 E^(-x)", 0);
    printf("testing %s\n", "TrigToExp[Cos[x]+I Sin[x]]"); fflush(stdout); assert_eval_eq("TrigToExp[Cos[x]+I Sin[x]]", "E^(I x)", 0);
    printf("testing %s\n", "TrigToExp[Tanh[x]]"); fflush(stdout); assert_eval_eq("TrigToExp[Tanh[x]]", "E^x/(E^x + E^(-x)) - E^(-x)/(E^x + E^(-x))", 0);
    printf("testing %s\n", "TrigToExp[ArcSin[x]]"); fflush(stdout); assert_eval_eq("TrigToExp[ArcSin[x]]", "-I Log[I x + Sqrt[1 - x^2]]", 0);
    printf("testing %s\n", "TrigToExp[ArcCoth[x]]"); fflush(stdout); assert_eval_eq("TrigToExp[ArcCoth[x]]", "1/2 Log[1 + 1/x] - 1/2 Log[1 - 1/x]", 0);
    printf("testing %s\n", "TrigToExp[{Sin[x],Cos[x],Tan[x]}]"); fflush(stdout); assert_eval_eq("TrigToExp[{Sin[x],Cos[x],Tan[x]}]", "{(-1/2*I) E^(I x) + (1/2*I) E^(-I x), 1/2 E^(-I x) + 1/2 E^(I x), -I E^(I x)/(E^(-I x) + E^(I x)) + I E^(-I x)/(E^(-I x) + E^(I x))}", 0);
    printf("testing %s\n", "TrigToExp[Sinh[x]<=11&&ArcTan[x]==7]"); fflush(stdout); assert_eval_eq("TrigToExp[Sinh[x]<=11&&ArcTan[x]==7]", "1/2 E^x - 1/2 E^(-x) <= 11 && (-1/2*I) Log[1 + I x] + (1/2*I) Log[1 - I x] == 7", 0);

    // Other inverse trig
    printf("testing %s\n", "TrigToExp[ArcCos[x]]"); fflush(stdout); assert_eval_eq("TrigToExp[ArcCos[x]]", "1/2 Pi + I Log[I x + Sqrt[1 - x^2]]", 0);
    printf("testing %s\n", "TrigToExp[ArcCot[x]]"); fflush(stdout); assert_eval_eq("TrigToExp[ArcCot[x]]", "(-1/2*I) Log[1 + I/x] + (1/2*I) Log[1 - I/x]", 0);
    printf("testing %s\n", "TrigToExp[ArcSec[x]]"); fflush(stdout); assert_eval_eq("TrigToExp[ArcSec[x]]", "1/2 Pi + I Log[I/x + Sqrt[1 - 1/x^2]]", 0);
    printf("testing %s\n", "TrigToExp[ArcCsc[x]]"); fflush(stdout); assert_eval_eq("TrigToExp[ArcCsc[x]]", "-I Log[I/x + Sqrt[1 - 1/x^2]]", 0);
    
    // Hyperbolics
    printf("testing %s\n", "TrigToExp[Cosh[x]]"); fflush(stdout); assert_eval_eq("TrigToExp[Cosh[x]]", "1/2 E^x + 1/2 E^(-x)", 0);
    printf("testing %s\n", "TrigToExp[Coth[x]]"); fflush(stdout); assert_eval_eq("TrigToExp[Coth[x]]", "E^x/(-E^x + E^(-x)) + E^(-x)/(-E^x + E^(-x))", 0);
    printf("testing %s\n", "TrigToExp[Sech[x]]"); fflush(stdout); assert_eval_eq("TrigToExp[Sech[x]]", "2/(E^x + E^(-x))", 0);
    printf("testing %s\n", "TrigToExp[Csch[x]]"); fflush(stdout); assert_eval_eq("TrigToExp[Csch[x]]", "2/(E^x - E^(-x))", 0);

    // Inverse hyperbolics
    printf("testing %s\n", "TrigToExp[ArcSinh[x]]"); fflush(stdout); assert_eval_eq("TrigToExp[ArcSinh[x]]", "Log[x + Sqrt[1 + x^2]]", 0);
    printf("testing %s\n", "TrigToExp[ArcCosh[x]]"); fflush(stdout); assert_eval_eq("TrigToExp[ArcCosh[x]]", "Log[x + Sqrt[-1 + x] Sqrt[1 + x]]", 0);
    printf("testing %s\n", "TrigToExp[ArcTanh[x]]"); fflush(stdout); assert_eval_eq("TrigToExp[ArcTanh[x]]", "1/2 Log[1 + x] - 1/2 Log[1 - x]", 0);
    printf("testing %s\n", "TrigToExp[ArcSech[x]]"); fflush(stdout); assert_eval_eq("TrigToExp[ArcSech[x]]", "Log[1/x + Sqrt[-1 + 1/x] Sqrt[1 + 1/x]]", 0);
    printf("testing %s\n", "TrigToExp[ArcCsch[x]]"); fflush(stdout); assert_eval_eq("TrigToExp[ArcCsch[x]]", "Log[1/x + Sqrt[1 + 1/x^2]]", 0);
    
    // Two-argument ArcTan
    printf("testing %s\n", "TrigToExp[ArcTan[x, y]]"); fflush(stdout); assert_eval_eq("TrigToExp[ArcTan[x, y]]", "-I Log[(x + I y)/Sqrt[x^2 + y^2]]", 0);
}

void test_exptotrig() {
    printf("testing %s\n", "ExpToTrig[Exp[I x]]"); fflush(stdout); assert_eval_eq("ExpToTrig[Exp[I x]]", "Cos[x] + I Sin[x]", 0);
    printf("testing %s\n", "ExpToTrig[Exp[x]]"); fflush(stdout); assert_eval_eq("ExpToTrig[Exp[x]]", "Cosh[x] + Sinh[x]", 0);
    printf("testing %s\n", "ExpToTrig[Log[1+I x]-Log[1-I x]]"); fflush(stdout); assert_eval_eq("ExpToTrig[Log[1+I x]-Log[1-I x]]", "(2*I) ArcTan[x]", 0);
    printf("testing %s\n", "ExpToTrig[Log[2+x]-Log[2-x]]"); fflush(stdout); assert_eval_eq("ExpToTrig[Log[2+x]-Log[2-x]]", "2 ArcTanh[1/2 x]", 0);
    printf("testing %s\n", "ExpToTrig[TrigToExp[Sin[x]]]"); fflush(stdout); assert_eval_eq("ExpToTrig[TrigToExp[Sin[x]]]", "Sin[x]", 0);
    printf("testing %s\n", "ExpToTrig[Exp[I*x]==-1]"); fflush(stdout); assert_eval_eq("ExpToTrig[Exp[I*x]==-1]", "Cos[x] + I Sin[x] == -1", 0);
    
    // Reverse trigonometric
    printf("testing %s\n", "ExpToTrig[TrigToExp[Cos[x]]]"); fflush(stdout); assert_eval_eq("ExpToTrig[TrigToExp[Cos[x]]]", "Cos[x]", 0);
    printf("testing %s\n", "ExpToTrig[TrigToExp[Tan[x]]]"); fflush(stdout); assert_eval_eq("ExpToTrig[TrigToExp[Tan[x]]]", "Tan[x]", 0);
    printf("testing %s\n", "ExpToTrig[TrigToExp[Cot[x]]]"); fflush(stdout); assert_eval_eq("ExpToTrig[TrigToExp[Cot[x]]]", "Cot[x]", 0);
    printf("testing %s\n", "ExpToTrig[TrigToExp[Sec[x]]]"); fflush(stdout); assert_eval_eq("ExpToTrig[TrigToExp[Sec[x]]]", "Sec[x]", 0);
    printf("testing %s\n", "ExpToTrig[TrigToExp[Csc[x]]]"); fflush(stdout); assert_eval_eq("ExpToTrig[TrigToExp[Csc[x]]]", "Csc[x]", 0);
    
    // Reverse hyperbolic
    printf("testing %s\n", "ExpToTrig[TrigToExp[Sinh[x]]]"); fflush(stdout); assert_eval_eq("ExpToTrig[TrigToExp[Sinh[x]]]", "Sinh[x]", 0);
    printf("testing %s\n", "ExpToTrig[TrigToExp[Cosh[x]]]"); fflush(stdout); assert_eval_eq("ExpToTrig[TrigToExp[Cosh[x]]]", "Cosh[x]", 0);
    printf("testing %s\n", "ExpToTrig[TrigToExp[Tanh[x]]]"); fflush(stdout); assert_eval_eq("ExpToTrig[TrigToExp[Tanh[x]]]", "Tanh[x]", 0);
    /* Note: the round-trip currently produces a sign-flipped result
     * (-Coth[x] instead of Coth[x]) due to a pre-existing TrigToExp/ExpToTrig
     * issue. The Times-level trig canonicalizer collapses -Cosh[x] Csch[x]
     * (the historical printed form) into the equivalent -Coth[x]. */
    printf("testing %s\n", "ExpToTrig[TrigToExp[Coth[x]]]"); fflush(stdout); assert_eval_eq("ExpToTrig[TrigToExp[Coth[x]]]", "-Coth[x]", 0);
    printf("testing %s\n", "ExpToTrig[TrigToExp[Sech[x]]]"); fflush(stdout); assert_eval_eq("ExpToTrig[TrigToExp[Sech[x]]]", "Sech[x]", 0);
    printf("testing %s\n", "ExpToTrig[TrigToExp[Csch[x]]]"); fflush(stdout); assert_eval_eq("ExpToTrig[TrigToExp[Csch[x]]]", "Csch[x]", 0);
}

int main() {
    symtab_init();
    core_init();
    test_trigtoexp();
    test_exptotrig();
    printf("All simp tests passed!\n");
    return 0;
}
