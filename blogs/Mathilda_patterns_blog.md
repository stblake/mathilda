# Patterns in Mathilda

In this short note I would like to explore pattern matching in Mathilda. One of the core features of the Mathematica programming language is its highly flexible and complex pattern matcher. 

Mathilda has the goal of supporting the core programming language of Mathematica, with the algebraic and calculus capabilities of (a long defunct) computer algebra system, Derive, all while having a minimal footprint [1]. 

Here is an extremely simple pattern in Mathilda, given $3 x^2$ we replace $x$ with $t$: 

```mathematica
In[1]:= 3 x^2 /. x -> t
Out[1]= 3 t^2
```

Let's consider a more complicated pattern where we replace all variables in a tri-variate polynomial by $t$. This makes use of the `Alternatives` primitive:

```mathematica
In[2]:= 3 x^3 + 2 y z /. (x | y | z) :> t
Out[2]= 2 t^2 + 3 t^3
```

Let's write a pattern matching rule that decreases the exponent of a variable (`Symbol`) by 1: 

```mathematica
In[3]:= A x^4 + B y^3 + C /. t_Symbol^n_Integer :> t^(n-1)
Out[3]= C + A x^3 + B y^2
```

Here's a pattern that uses `PatternTest` to check if an integer is even: 

```mathematica
In[4]:= Range[10] /. n_Integer?EvenQ :> n/2
Out[4]= {1, 1, 3, 2, 5, 3, 7, 4, 9, 5}
```

And here we mix ReplaceRepeated (`//.`) and Conditions on patterns to implement a check for the Collatz Conjecture. This check requires implementing two rules based on the parity of the number: 

```mathematica
(* In Mathilda trailing slashes (\) are used used for repl line continuation, a la Python. *)
In[5]:= collatzRules = {\
  {pre___, n_Integer} /; EvenQ[n] && n > 1 :> {pre, n, n / 2},\
  {pre___, n_Integer} /; OddQ[n] && n > 1 :> {pre, n, 3 n + 1}\
};

In[6]:= {27} //. collatzRules
Out[6]= {27, 82, 41, 124, 62, 31, 94, 47, 142, 71, 214, 107, 322, 161, 484, 242, 121, 364, 182, 91, 274, 137, 412, 206, 103, 310, 155, 466, 233, 700, 350, 175, 526, 263, 790, 395, 1186, 593, 1780, 890, 445, 1336, 668, 334, 167, 502, 251, 754, 377, 1132, 566, 283, 850, 425, 1276, 638, 319, 958, 479, 1438, 719, 2158, 1079, 3238, 1619, 4858, 2429, 7288, 3644, 1822, 911, 2734, 1367, 4102, 2051, 6154, 3077, 9232, 4616, 2308, 1154, 577, 1732, 866, 433, 1300, 650, 325, 976, 488, 244, 122, 61, 184, 92, 46, 23, 70, 35, 106, 53, 160, 80, 40, 20, 10, 5, 16, 8, 4, 2, 1}

In[7]:= Max[%]
Out[7]= 9232

In[8]:= Length[%%]
Out[8]= 112
```
We see that starting with 27, the sequence increases up to 9232 before finally settling to 1 after 112 iterations. 

Next we will highlight the use of multiple `BlankNullSequence` patterns and the backtracking feature of the pattern matcher by implementing a rule that removes duplicate entries in a list: 

```mathematica
In[9]:= deldups = {a___, x_, b___, x_, c___} :> {a, x, b, c};

In[10]:= Mod[Range[30]^2,17]
Out[10]= {1, 4, 9, 16, 8, 2, 15, 13, 13, 15, 2, 8, 16, 9, 4, 1, 0, 1, 4, 9, 16, 8, 2, 15, 13, 13, 15, 2, 8, 16}

In[11]:= % //. deldups
Out[11]= {1, 4, 9, 16, 8, 2, 15, 13, 0}
```

Sorting a list usually requires loops and temporary variables. With patterns, sorting is just a rule: If you see two adjacent elements out of order, swap them:

```mathematica
In[12]:= lst = {2, 3, 6, 8, 10, 6, 7, 6, 9, 10, 3, 4, 5, 3, 6, 4, 6, 3, 6, 1};

In[13]:= lst //. {pre___, a_, b_, post___} /; a > b :> {pre, b, a, post}
Out[13]= {1, 2, 3, 3, 3, 3, 4, 4, 5, 6, 6, 6, 6, 6, 6, 7, 8, 9, 10, 10}
```

Similarly, we can make a one line implementation of the Sieve of Eratosthenes: 

```mathematica
In[14]:= Range[2, 100] //. {pre___, x_, mid___, y_, post___} /; Mod[y, x] == 0 :> {pre, x, mid, post}
Out[14]= {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97}
```

It's easy to confirm they are all prime numbers: 
```
In[15]:= PrimeQ /@ %
Out[15]= {True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True}
```

Below we see a more complex pattern that implements a run-length encoding. In addition to `BlankNullSequence`, `Blank`, and `Condition`, this pattern makes use of `Longest` and `Repeated`: 

```mathematica
In[16]:= rleRule = {pre___, r:Longest[seq_ ..], post___} /; Length[{r}]>1 :> {pre, {seq, Length[{r}]}, post};

In[17]:= data = {"A", "A", "A", "B", "B", "C", "A", "A"};

In[18]:= data //. rleRule 
Out[18]= {{"A", 3}, {"B", 2}, "C", {"A", 2}}
```

Patterns are not just limited to transformation rules, here's an example where we count the number of primes less than 1000 that are of the form 4 n + 3: 

```mathematica
In[19]:= Count[Range[1000], p_?PrimeQ /; Mod[p,4] == 3]
Out[19]= 87
```

Next we'll look at an application of `Flat` pattern matching. As multiplication is associative, it has the attribute `Flat`:

```mathematica
In[20]:= MemberQ[Attributes[Times], Flat] 
Out[20]= True
```

We can implement a rule for expanding out the logarithm of products: 

```mathematica
In[21]:= logExpand = Log[x_ * y_] :> Log[x] + Log[y];

In[22]:= Log[a * b * c * d] /. logExpand
Out[22]= Log[a] + Log[b c d]
```

This works because the pattern matcher matches `x_ y_` to `Times[a, Times[b, c, d]`. We can apply this rule two more times to fully expand out the expression: 

```mathematica
In[23]:= % /. logExpand
Out[23]= Log[a] + Log[b] + Log[c d]

In[24]:= % /. logExpand
Out[24]= Log[a] + Log[b] + Log[c] + Log[d]
```

Alternatively, we could have done this with `ReplaceRepeated` (`//.`): 

```mathematica
In[25]:= Log[a * b * c * d] //. logExpand
Out[25]= Log[a] + Log[b] + Log[c] + Log[d]
```

The pattern matcher supports backtracking, here's an example where we find all positive integers `a, b, c`, such that `a + b c == 71`:

```mathematica
In[26]:= ReplaceList[Range[71], {a___,x_,b___,y_,c___,z_,d___} /; x + y z === 71 :> {x,y,z}]
Out[26]= {{1, 2, 35}, {1, 5, 14}, {1, 7, 10}, {2, 3, 23}, {3, 4, 17}, {5, 6, 11}}

In[27]:= Print[HoldForm[ #1 + #2 #3 == 71 ] // InputForm] & @@@ %;
1 + 2 35 == 71
1 + 5 14 == 71
1 + 7 10 == 71
2 + 3 23 == 71
3 + 4 17 == 71
5 + 6 11 == 71
```

We can also use pattern constraints `Longest` and `Shortest` to find the smallest and largest difference between `x` and `y`:

```mathematica
In[28]:= Range[71] /. {a___,x_,Shortest[b___],y_,c___,z_,d___} /; x + y z === 71 :> {x,y,z}
Out[28]= {1, 2, 35}

In[29]:= Range[71] /. {a___,x_,Longest[b___],y_,c___,z_,d___} /; x + y z === 71 :> {x,y,z}
Out[29]= {1, 7, 10}
```

[1] https://en.wikipedia.org/wiki/Derive_(computer_algebra_system)

