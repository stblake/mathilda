
# Mathilda — The first AI Agent Generated Computer Algebra System

Sam Blake, April 2026

## Building Mathilda

Since I was an undergraduate I have loved the Mathematica computer algebra system. It’s a brilliant language and great system for doing both mathematics and general research by computer. 

Given the rapidly evolving landscape of LLM-based coding agents, I wondered if now is the right time to try an AI Agent build of a small computer algebra system. My first thought was to recreate the most basic functionality of Mathematica using Gemini 3 pro CLI. This would require some significant software development, including: 

1. The necessary data structures to hold symbolic expressions in a low level language like C. 
2. Parsing a complex language like Mathematica into the directed graph data structure of a symbolic expression. 
3. Constructing a hash table for storing the symbolic expressions. 
4. Building a pattern matcher capable of matching Blank, BlankSequence, PatternTest, Condition etc. with backtracking. (Including matching Orderless, Flat expressions etc.)
5. Building a recreation of the infinite evaluation semantics of Mathematica, including Attributes (Flat, Listable, Orderless etc,) unevaluated expressions (Hold, HoldPattern etc,) evaluation based on DownValues, UpValues, OwnValues etc. 
6. Building the minimal algebraic capabilities expected of a tiny computer algebra system. 

I really wasn’t sure an AI Agent was up to it in March, 2026. I couldn’t have been more wrong. Within 3 days of starting this project with minimal oversight (I was away on a golf trip while the initial build took place) and without writing a single line of code, we could evaluate expressions like: 

```
$ ./Mathilda 

Mathilda - A tiny, LLM-generated, Mathematica-like computer algebra system.

This program is free, open source software and comes with ABSOLUTELY NO WARRANTY.

End a line with '\' to enter a multiline expression. Press Return to evaluate.
Exit by evaluating Quit[] or CONTROL-C.

In[1]:= lst = {1,3,4,2,3,4,5,4,2,4,6,2,1}
Out[1]= {1, 3, 4, 2, 3, 4, 5, 4, 2, 4, 6, 2, 1}

In[2]:= lst //. {a___, x_, b___, x_, c___} -> {a, x, b, c}
Out[2]= {1, 3, 4, 2, 5, 6}

In[3]:= UnsortedUnion[x_] := Module[{f}, f[y_] := (f[y] = Sequence[]; y); f /@ x]
Out[3]= Null

In[4]:= lst // UnsortedUnion
Out[4]= {1, 3, 4, 2, 5, 6}

In[5]:= Expand[(x + 2 y)^3 (3 x + 4 y)^2]
Out[5]= 384 x y^4 + 9 x^5 + 456 x^2 y^3 + 268 x^3 y^2 + 78 x^4 y + 128 y^5

In[6]:= % // FactorSquareFree
Out[6]= (x + 2 y)^3 (3 x + 4 y)^2

In[7]:= {x, y, z} . {{a, b, c}, {d, e, f}, {g, h, i}} . {p, q, r}
Out[7]= p (a x + d y + g z) + q (b x + e y + h z) + r (c x + f y + i z)
```

To me, this was a very impressive result. Gemini 3 pro CLI has built the basics of a complex programming language from scratch in a low level language, C. Thinking about it more, none of the individual parts of the build are especially challenging for a mid/senior dev (with a vast general knowledge of computer algebra design and a detailed knowledge of the Mathematica programming language.) All of the data structures and algorithms required are well known and the Mathematica language is well documented on the internet, but the ability to put all these pieces together so quickly was staggering! 

I find myself wondering what Mathilda would look like with several AI Agent years of development. (At present it has had significantly less than a week of continuous agentic development.) I highly doubt it will be implementing the algebraic case of the Risch-Trager algorithm anytime soon, but I would be happily surprised if it did. 


## AI Agent Code Development Takeaways 

Gemini 3 pro CLI produces significantly better code than Gemini 3 Flash CLI. Due to Flash induced memory leaks, crashes, and bugs; the development time for pro is significantly faster and more cost effective. 

Chunking up tasks into manageable, discrete, well defined pieces vastly improves development time and code quality. 

Strictly enforcing the creation and running of automated unit tests (via cmake) and checking for memory leaks (via valgrind) after each change drastically improves code quality. 

Automating commits for discrete improvements via Gemini establishes reliable rollback points in the event of feature regressions.


## Mathilda v Mathematica

Mathilda does not remotely compare to Mathematica, neither in performance, functionality, documentation, battle-hardened code etc. At best Mathilda could be jovially thought of as a gateway drug to Mathematica. Given its open source licence, Mathilda is a good learning tool for students and researchers who want to look under the hood and learn how a bare-bones computer algebra system is put together in a low level language like C. 


## The Name - Mathilda

A quick note on the name, Mathilda. As Mathilda is designed to be small (tiny by modern standards,) with a compiled binary of just a few hundred kB. I thought a name that conveyed its goal of being small (pico being the prefix for 1e-12) and easily extensible by both human and AI developers made sense. It’s also in homage to the tiny computer algebra systems of the 1980s by David Stoutemyer and Albert Rich PICOMATH, MuMath, and the brilliant, yet tiny CAS, Derive. 

