#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>

void test_distribute() {
    // Apply the distributive law:
    assert_eval_eq("Distribute[(a+b).(x+y+z)]", "Dot[a, x] + Dot[b, x] + Dot[a, y] + Dot[b, y] + Dot[a, z] + Dot[b, z]", 0);
    
    // Distribute f over Plus:
    assert_eval_eq("Distribute[f[a+b,c+d+e]]", "f[a, c] + f[b, c] + f[a, d] + f[b, d] + f[a, e] + f[b, e]", 0);
    
    // Distribute f over g:
    assert_eval_eq("Distribute[f[g[a,b],g[c,d,e]],g]", "g[f[a, c], f[a, d], f[a, e], f[b, c], f[b, d], f[b, e]]", 0);
    
    // By default, distribute over Plus:
    assert_eval_eq("Distribute[(a+b+c)(u+v),Plus]", "a u + b u + c u + a v + b v + c v", 0);
    assert_eval_eq("Distribute[(a+b+c)(u+v)]", "a u + b u + c u + a v + b v + c v", 0);
    
    // Distribute only a product of sums:
    assert_eval_eq("Distribute[(a+b+c)(u+v),Plus,Times]", "a u + b u + c u + a v + b v + c v", 0);
    assert_eval_eq("Distribute[(a+b+c)^(u+v),Plus,Times]", "(a + b + c)^(u + v)", 0);
    
    // Distribute f over g to give fp and gp:
    assert_eval_eq("Distribute[f[g[a,b],g[c,d,e]],g,f,gp,fp]", "gp[fp[a, c], fp[a, d], fp[a, e], fp[b, c], fp[b, d], fp[b, e]]", 0);
    assert_eval_eq("Distribute[(a+b+c)(u+v),Plus,Times,plus,times]", "plus[times[a, u], times[a, v], times[b, u], times[b, v], times[c, u], times[c, v]]", 0);
    
    // Expand symbolic matrix and vector expressions:
    assert_eval_eq("Distribute[(m1+m2).(m3+m4)]", "Dot[m1, m3] + Dot[m2, m3] + Dot[m1, m4] + Dot[m2, m4]", 0);
    assert_eval_eq("Distribute[Cross[v1+v2, v3+v4]]", "Cross[v1, v3] + Cross[v2, v3] + Cross[v1, v4] + Cross[v2, v4]", 0);
    
    // Find the list of all possible combinations of elements:
    assert_eval_eq("Distribute[{{a, b}, {x,y,z},{s,t}}, List]", "{{a, x, s}, {a, x, t}, {a, y, s}, {a, y, t}, {a, z, s}, {a, z, t}, {b, x, s}, {b, x, t}, {b, y, s}, {b, y, t}, {b, z, s}, {b, z, t}}", 0);
    
    // Generate a powerset:
    assert_eval_eq("Distribute[{{{},{a}}, {{},{b}}, {{},{c}}}, List, List, List, Join]", "{{}, {c}, {b}, {b, c}, {a}, {a, c}, {a, b}, {a, b, c}}", 0);

    // Find intermediate terms from a direct application of the distributive law:
    assert_eval_eq("Distribute[(-1 + x) (1 + x) (1 - x + x^2) (1 + x + x^2), Plus, Times, List, Times]", 
        "{-1, x, -x^2, -x, x^2, -x^3, -x^2, x^3, -x^4, -x, x^2, -x^3, -x^2, x^3, -x^4, -x^3, x^4, -x^5, x, -x^2, x^3, x^2, -x^3, x^4, x^3, -x^4, x^5, x^2, -x^3, x^4, x^3, -x^4, x^5, x^4, -x^5, x^6}", 0);

    // For pure products, Distribute gives the same results as Expand:
    assert_eval_eq("Distribute[(-1 + x) (1 + x) (1 - x + x^2 - x^3 + x^4) (1 + x + x^2 + x^3 + x^4)]", "-1 + x^10", 0);
}

void test_inner() {
    assert_eval_eq("Inner[f,{a,b},{x,y},g]", "g[f[a, x], f[b, y]]", 0);
    assert_eval_eq("Inner[f,{{a,b},{c,d}},{x,y},g]", "{g[f[a, x], f[b, y]], g[f[c, x], f[d, y]]}", 0);
    assert_eval_eq("Inner[Times,{a,b},{x,y},Plus]", "a x + b y", 0);
    assert_eval_eq("Inner[Power,{a,b,c},{x,y,z},Times]", "a^x b^y c^z", 0);
    assert_eval_eq("Inner[f,{{a,b},{c,d}},{{u,v},{w,x}},g]", "{{g[f[a, u], f[b, w]], g[f[a, v], f[b, x]]}, {g[f[c, u], f[d, w]], g[f[c, v], f[d, x]]}}", 0);
    assert_eval_eq("Inner[f,{{a,b},{c,d}},{s,t},g]", "{g[f[a, s], f[b, t]], g[f[c, s], f[d, t]]}", 0);
    assert_eval_eq("Inner[f,{x,y},{{a,b},{c,d}},g]", "{g[f[x, a], f[y, c]], g[f[x, b], f[y, d]]}", 0);
    assert_eval_eq("Inner[f,{{a,b},{c,d}},{{u,v},{w,x}},g,1]", "{{g[f[a, u], f[c, w]], g[f[a, v], f[c, x]]}, {g[f[b, u], f[d, w]], g[f[b, v], f[d, x]]}}", 0);
    assert_eval_eq("Inner[Times,{a1, a2, a3},{b1, b2, b3},Plus]", "a1 b1 + a2 b2 + a3 b3", 0);
}

void test_outer() {
    assert_eval_eq("Outer[f,{a,b},{x,y,z}]", "{{f[a, x], f[a, y], f[a, z]}, {f[b, x], f[b, y], f[b, z]}}", 0);
    assert_eval_eq("Outer[Times,{1,2,3,4},{a,b,c}]", "{{a, b, c}, {2 a, 2 b, 2 c}, {3 a, 3 b, 3 c}, {4 a, 4 b, 4 c}}", 0);
    assert_eval_eq("Outer[Times,{{1,2},{3,4}},{{a,b},{c,d}}]", "{{{{a, b}, {c, d}}, {{2 a, 2 b}, {2 c, 2 d}}}, {{{3 a, 3 b}, {3 c, 3 d}}, {{4 a, 4 b}, {4 c, 4 d}}}}", 0);
    assert_eval_eq("Outer[f,{a,b},{x,y,z},{u,v}]", "{{{f[a, x, u], f[a, x, v]}, {f[a, y, u], f[a, y, v]}, {f[a, z, u], f[a, z, v]}}, {{f[b, x, u], f[b, x, v]}, {f[b, y, u], f[b, y, v]}, {f[b, z, u], f[b, z, v]}}}", 0);
    assert_eval_eq("Outer[f,{{1,2},{3,4}},{{a,b},{c,d}},1]", "{{f[{1, 2}, {a, b}], f[{1, 2}, {c, d}]}, {f[{3, 4}, {a, b}], f[{3, 4}, {c, d}]}}", 0);
    assert_eval_eq("Outer[Times,{{1,2},{3,4}},{{a,b,c},{d,e}}]", "{{{{a, b, c}, {d, e}}, {{2 a, 2 b, 2 c}, {2 d, 2 e}}}, {{{3 a, 3 b, 3 c}, {3 d, 3 e}}, {{4 a, 4 b, 4 c}, {4 d, 4 e}}}}", 0);
    assert_eval_eq("Outer[g,f[a,b],f[x,y,z]]", "f[f[g[a, x], g[a, y], g[a, z]], f[g[b, x], g[b, y], g[b, z]]]", 0);
    assert_eval_eq("Outer[f]", "f[]", 0);
    assert_eval_eq("Outer[List,{a,b,c},{x,y}]", "{{{a, x}, {a, y}}, {{b, x}, {b, y}}, {{c, x}, {c, y}}}", 0);
}

void test_tuples() {
    assert_eval_eq("Tuples[{0,1},3]", "{{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}}", 0);
    assert_eval_eq("Tuples[{1,0},3]", "{{1, 1, 1}, {1, 1, 0}, {1, 0, 1}, {1, 0, 0}, {0, 1, 1}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}}", 0);
    assert_eval_eq("Tuples[{{a,b},{1,2,3,4},{x}}]", "{{a, 1, x}, {a, 2, x}, {a, 3, x}, {a, 4, x}, {b, 1, x}, {b, 2, x}, {b, 3, x}, {b, 4, x}}", 0);
    assert_eval_eq("Tuples[{a,b},2]", "{{a, a}, {a, b}, {b, a}, {b, b}}", 0);
    assert_eval_eq("Tuples[{a,a,b},2]", "{{a, a}, {a, a}, {a, b}, {a, a}, {a, a}, {a, b}, {b, a}, {b, a}, {b, b}}", 0);
    assert_eval_eq("Tuples[{a,b},{2,2}]", "{{{a, a}, {a, a}}, {{a, a}, {a, b}}, {{a, a}, {b, a}}, {{a, a}, {b, b}}, {{a, b}, {a, a}}, {{a, b}, {a, b}}, {{a, b}, {b, a}}, {{a, b}, {b, b}}, {{b, a}, {a, a}}, {{b, a}, {a, b}}, {{b, a}, {b, a}}, {{b, a}, {b, b}}, {{b, b}, {a, a}}, {{b, b}, {a, b}}, {{b, b}, {b, a}}, {{b, b}, {b, b}}}", 0);
    assert_eval_eq("Tuples[f[x,y,z],2]", "{f[x, x], f[x, y], f[x, z], f[y, x], f[y, y], f[y, z], f[z, x], f[z, y], f[z, z]}", 0);
    assert_eval_eq("Tuples[{{a,b},{x,y,z}}]", "{{a, x}, {a, y}, {a, z}, {b, x}, {b, y}, {b, z}}", 0);
    assert_eval_eq("Tuples[{{a,b},{x,y,z}},2]", "{{{a, b}, {a, b}}, {{a, b}, {x, y, z}}, {{x, y, z}, {a, b}}, {{x, y, z}, {x, y, z}}}", 0);
}

void test_permutations() {
    assert_eval_eq("Permutations[{a,b,c}]", "{{a, b, c}, {a, c, b}, {b, a, c}, {b, c, a}, {c, a, b}, {c, b, a}}", 0);
    assert_eval_eq("Permutations[{a,b,c,d},{3}]", "{{a, b, c}, {a, b, d}, {a, c, b}, {a, c, d}, {a, d, b}, {a, d, c}, {b, a, c}, {b, a, d}, {b, c, a}, {b, c, d}, {b, d, a}, {b, d, c}, {c, a, b}, {c, a, d}, {c, b, a}, {c, b, d}, {c, d, a}, {c, d, b}, {d, a, b}, {d, a, c}, {d, b, a}, {d, b, c}, {d, c, a}, {d, c, b}}", 0);
    assert_eval_eq("Permutations[{a,a,b}]", "{{a, a, b}, {a, b, a}, {b, a, a}}", 0);
    assert_eval_eq("Permutations[{x,x^2,x+1}]", "{{x, x^2, 1 + x}, {x, 1 + x, x^2}, {x^2, x, 1 + x}, {x^2, 1 + x, x}, {1 + x, x, x^2}, {1 + x, x^2, x}}", 0);
    assert_eval_eq("Permutations[Range[3],All]", "{{}, {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 1}, {2, 3}, {3, 1}, {3, 2}, {1, 2, 3}, {1, 3, 2}, {2, 1, 3}, {2, 3, 1}, {3, 1, 2}, {3, 2, 1}}", 0);
    assert_eval_eq("Permutations[Range[4],{4,0,-2}]", "{{1, 2, 3, 4}, {1, 2, 4, 3}, {1, 3, 2, 4}, {1, 3, 4, 2}, {1, 4, 2, 3}, {1, 4, 3, 2}, {2, 1, 3, 4}, {2, 1, 4, 3}, {2, 3, 1, 4}, {2, 3, 4, 1}, {2, 4, 1, 3}, {2, 4, 3, 1}, {3, 1, 2, 4}, {3, 1, 4, 2}, {3, 2, 1, 4}, {3, 2, 4, 1}, {3, 4, 1, 2}, {3, 4, 2, 1}, {4, 1, 2, 3}, {4, 1, 3, 2}, {4, 2, 1, 3}, {4, 2, 3, 1}, {4, 3, 1, 2}, {4, 3, 2, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 1}, {2, 3}, {2, 4}, {3, 1}, {3, 2}, {3, 4}, {4, 1}, {4, 2}, {4, 3}, {}}", 0);
    assert_eval_eq("Permutations[f[a,b,c]]", "{f[a, b, c], f[a, c, b], f[b, a, c], f[b, c, a], f[c, a, b], f[c, b, a]}", 0);
}

int main() {
    symtab_init();
    core_init();
    test_distribute();
    test_inner();
    test_outer();
    test_tuples();
    test_permutations();
    printf("All distribute tests passed!\\n");
    return 0;
}
