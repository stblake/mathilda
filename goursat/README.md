# Goursat (1887): Pseudo-Elliptic Integrals — Review and Implementation

**Reference.** É. Goursat, *Note sur quelques intégrales pseudo-elliptiques*, Bulletin de la S.M.F., t. 15 (1887), pp. 106–120 (séance du 16 mars 1887).

This package implements the algorithmic content of Goursat's note in Mathematica. The two files are

- `GoursatPseudoElliptic.wl` — the package proper.
- `GoursatExamples.wl` — six worked examples (positive, negative, finite-fixed-point, symbolic, cubic).

Load the package with `Get["GoursatPseudoElliptic.wl"]`.

---

## 1. Review of the paper

A *pseudo-elliptic integral* is an integral whose integrand looks elliptic — that is, of the form

$$\int \mathcal{F}(t)\,\frac{dt}{\sqrt{R(t)}}$$

with $\mathcal{F}$ rational and $R$ a polynomial of degree 3 or 4 with simple roots — but which happens to evaluate in elementary terms. Goursat's note characterises a wide class of such integrals via the *Möbius involutions that permute the roots of $R$*.

### The forward construction (§1)

Goursat's first observation is constructive. Take any elementary integrand
$$f(x)\,\frac{dx}{\sqrt{Ax^2+Bx+C}}$$
and apply a quadratic-rational change of variable $x=\varphi(t)=(at^2+bt+c)/(a't^2+b't+c')$. The result is a new integrand of the form $\mathcal{F}(t)\,dt/\sqrt{R(t)}$ with $R$ a quartic. The key algebraic fact, established by direct computation in §1, is that the equation $\varphi(t)=\lambda$ has two roots $t_1,t_2$ tied by a *fixed* involutive linear relation
$$L\,t_1 t_2 - N(t_1+t_2) + M = 0, \qquad L=ab'\!-\!ba',\ M=bc'\!-\!cb',\ N=ca'\!-\!ac'.$$
Equivalently, the involution $S\colon t\mapsto (Nt-M)/(Lt-N)$ permutes $t_1,t_2$, and one verifies
$$\mathcal{F}(S(t)) = -\mathcal{F}(t).$$
So *every integrand obtained from an elementary one by this kind of substitution is anti-invariant under some Möbius involution* permuting the roots of $R$.

### Theorem 1 — the converse (§2)

Goursat then proves the converse, which is the heart of the matter:

> **Theorem 1.** If $S$ is a non-trivial Möbius involution that permutes the four roots of a quartic $R(t)$, and $\mathcal{F}$ is rational with $\mathcal{F}(t)+\mathcal{F}(S(t))=0$, then $\int \mathcal{F}\,dt/\sqrt{R}$ is elementary.

The proof is a beautiful two-step reduction. Let $\alpha,\beta$ be the fixed points of $S$. Substitute
$$u = \frac{t-\alpha}{t-\beta}.$$
Under this Möbius change, $S$ is conjugated to the simple reflection $u\mapsto -u$ (verify: $u_2 = (t_2-\alpha)/(t_2-\beta) = -u_1$ exactly when $S(t_1)=t_2$). Three things happen simultaneously:

1. $\mathcal{F}(t(u))$ is *odd* in $u$ (because anti-invariance under $S$ becomes anti-invariance under $u\mapsto -u$).
2. $R(t(u))\cdot(u-1)^4$ is *even* in $u$ (because $S$ permutes the roots of $R$ in $\pm$-pairs in $u$-coordinates).
3. The Jacobian factor $dt/du = (\alpha-\beta)/(u-1)^2$ collaborates with the $(u-1)^4$ inside the radical so that the differential
   $$\frac{\mathcal{F}(t)\,dt}{\sqrt{R(t)}} \;=\; \frac{(\alpha-\beta)\,\mathcal{F}(t(u))}{\sqrt{R(t(u))(u-1)^4}}\,du$$
   has *odd* numerator and *even* radicand in $u$.

A second substitution $x=u^2$, $u\,du = dx/2$, then turns this into
$$\frac{1}{2}\,\frac{G(x)}{\sqrt{H(x)}}\,dx$$
where $H(x)$ is at most quadratic in $x$ — manifestly elementary.

The special case where one fixed point is at $\infty$ collapses: the substitution becomes $u=t-\alpha$ directly, and the same argument runs through.

### Theorem 2 — the four-root sum criterion (§3)

For a generic quartic, there are exactly **three** non-trivial Möbius involutions permuting the four roots, corresponding to the three pairings of $\{a,b,c,d\}$ into two pairs. Goursat writes one down explicitly:
$$S_{(a,b)|(c,d)}(t) = \frac{(ab-cd)\,t + (a+b)cd - (c+d)ab}{[(a+b)-(c+d)]\,t - (ab-cd)}.$$
Together with the identity, $\{I,S_1,S_2,S_3\}$ form a Klein four-group $V_4 \subset \mathrm{PGL}_2$ (one checks $S_1S_2=S_3$, etc.).

> **Theorem 2.** If $\mathcal{F}$ is rational and $\mathcal{F}+\mathcal{F}(S_1)+\mathcal{F}(S_2)+\mathcal{F}(S_3)=0$, then $\int \mathcal{F}\,dt/\sqrt R$ is elementary.

The proof uses character theory of $V_4$. The condition says the projection of $\mathcal{F}$ onto the *trivial* character of $V_4$ vanishes, so $\mathcal{F}$ decomposes into projections onto the three non-trivial characters. Each non-trivial character is anti-invariant under at least one $S_i$, so each projection is Theorem-1-pseudo-elliptic. Adding the three antiderivatives gives the answer.

Goursat also notes (§3, end) that the sum $\mathcal{F}+\sum_i\mathcal{F}(S_i)$ is a symmetric function of the roots $a,b,c,d$, so the criterion can in principle be tested using the elementary symmetric functions (i.e. the coefficients of $R$) without ever computing the roots — a useful remark for symbolic verification.

### Sections 4–6 — special $R$ with extra symmetry

When $R$ has a special form, the group of Möbius transformations permuting its roots can be larger than $V_4$:

- **§4** — period-3 substitutions, e.g. $R(t) = t^3-1$, with cyclic action on $\{1,\omega,\omega^2\}$ fixing $\infty$. The hypothesis $\mathcal{F}(S(t)) = \omega\,\mathcal{F}(t)$ (with $S$ of order 3) gives a new family. The classical *Legendre / Clausen integral* $\int t\,dt/((t^3+8)\sqrt{t^3-1})$ falls in this case.

- **§5** — full tetrahedral group $A_4$ (12 elements) for $R=t^3-1$. The criterion involves two character-theoretic sums.

- **§6** — octahedral group for $R=t^4-1$ or, equivalently up to a Möbius change, $R=t(t^2+1)$. These give 8 extra substitutions, organised by a $V_4$ subgroup plus four order-4 elements.

- **§7** — higher-genus extensions (sextic $R$), connecting to Jacobi's reduction (Crelle 8) and a result of Poincaré on integrals of the first kind.

### A small correction

The blue-box "proof of Theorem 2" in the manuscript states that each projected component $\mathcal{F}_i$ is anti-invariant under $S_i$ with *matching index*. A direct character calculation shows this labelling is off:

| Component | Defining signs $(I,S_1,S_2,S_3)$ | Anti-inv. under | Invariant under |
|-----------|-----|-----|-----|
| $\mathcal{F}_1 = (\mathcal{F}-\mathcal{F}(S_1)-\mathcal{F}(S_2)+\mathcal{F}(S_3))/4$ | $(+,-,-,+)$ | $S_1, S_2$ | $S_3$ |
| $\mathcal{F}_2 = (\mathcal{F}-\mathcal{F}(S_1)+\mathcal{F}(S_2)-\mathcal{F}(S_3))/4$ | $(+,-,+,-)$ | $S_1, S_3$ | $S_2$ |
| $\mathcal{F}_3 = (\mathcal{F}+\mathcal{F}(S_1)-\mathcal{F}(S_2)-\mathcal{F}(S_3))/4$ | $(+,+,-,-)$ | $S_2, S_3$ | $S_1$ |

So $\mathcal{F}_2$ is *invariant* under $S_2$ (not anti-invariant). The conclusion of the theorem is unchanged — every $\mathcal{F}_i$ *is* anti-invariant under at least one involution — but the implementation must pick a valid involution for each component. The package uses $S_1$ for $\mathcal{F}_1$ and $\mathcal{F}_2$, and $S_2$ for $\mathcal{F}_3$.

---

## 2. The implementation

### Public functions

| Function | Purpose |
|-----------|---------|
| `PairSwapInvolution[a, b, c, d, t]` | Goursat's explicit formula for the involution swapping $\{a,b\}$ with $\{c,d\}$. Handles `Infinity` as one of the four points. |
| `QuarticRoots[R, t]` | Roots of a cubic or quartic; pads with `Infinity` if cubic. |
| `GoursatInvolutions[R, t]` | The three involutions $S_1, S_2, S_3$ for a given $R$. |
| `InvolutionFixedPoints[S, t]` | Fixed points of a Möbius involution; one may be `Infinity`. |
| `AntiInvariantQ[F, S, t]` | Theorem 1 hypothesis: $F + F(S) \stackrel{?}{=} 0$. |
| `PseudoEllipticQ[F, R, t]` | Theorem 2 hypothesis: $F + F(S_1) + F(S_2) + F(S_3) \stackrel{?}{=} 0$. |
| `GoursatReduceTheorem1[F, R, S, t]` | Performs the two-step substitution reduction and returns an antiderivative in $t$. |
| `GoursatPseudoElliptic[F, R, t]` | Top-level: tests the Theorem 2 hypothesis, decomposes $F$ into character components, and applies Theorem 1 to each. |

### Algorithm flow

```
GoursatPseudoElliptic[F, R, t]
   ├── GoursatInvolutions[R, t]       (three Möbius involutions)
   ├── verify F + F(S1) + F(S2) + F(S3) == 0
   ├── decompose F = F1 + F2 + F3 via V4 characters
   └── for each Fi, GoursatReduceTheorem1[Fi, R, Si, t]:
          ├── InvolutionFixedPoints[Si, t]
          ├── substitute u = (t-α)/(t-β)
          ├── verify integrand is odd in u over even-in-u radicand
          ├── substitute x = u^2  →  G(x)/(2 √H(x)) dx
          ├── Integrate[…, x]
          └── back-substitute x ← ((t-α)/(t-β))^2
```

### Worked example traced by hand

Take $R(t)=(t^2-1)(t^2-4)$, $F(t)=t$.

The roots are $\{\pm 1, \pm 2\}$ and the three involutions come out to (modulo the order Mathematica chooses)
$$S_1(t) = -t, \quad S_2(t) = 2/t, \quad S_3(t) = -2/t.$$
Then $F + F(S_1) + F(S_2) + F(S_3) = t + (-t) + 2/t + (-2/t) = 0$. ✓

Decomposing,
$$\mathcal{F}_1 = 0, \quad \mathcal{F}_2 = \frac{t}{2}-\frac{1}{t} = \frac{t^2-2}{2t}, \quad \mathcal{F}_3 = \frac{t}{2}+\frac{1}{t} = \frac{t^2+2}{2t}.$$

For $\mathcal{F}_3$ with $S_2 = 2/t$ (fixed points $\pm\sqrt 2$): the substitutions yield, after some algebra,
$$\int \frac{t^2+2}{2t}\frac{dt}{\sqrt{(t^2-1)(t^2-4)}} \;=\; -8\!\int\!\frac{dx}{(1-x)\sqrt{-2(x^2-34x+1)}}\bigg|_{x=((t-\sqrt 2)/(t+\sqrt 2))^2}$$
which is an arcsinh / log expression in elementary functions.

For $\mathcal{F}_2$ we use $S_1=-t$ (fixed points $\{0,\infty\}$): the substitution is just $u=t$, then $x=t^2$, and the integral collapses to $\int (x-2)/(4x\sqrt{x^2-5x+4})\,dx$.

Adding the two contributions and simplifying produces the closed form
$$\int \frac{t\,dt}{\sqrt{(t^2-1)(t^2-4)}} \;=\; \tfrac{1}{2}\log\!\big(2t^2 - 5 + 2\sqrt{(t^2-1)(t^2-4)}\big) + C,$$
which one verifies by direct differentiation.

### Limitations and caveats

1. **Only $V_4$.** The package implements Theorems 1 and 2. Sections 4–6 of the paper, including the Legendre/Clausen integral $\int t\,dt/((t^3+8)\sqrt{t^3-1})$, require the period-3 / tetrahedral / octahedral extensions. Adding them would mean implementing character projections for the larger groups $A_4$ and $S_4$ acting on the roots.

2. **Roots of $R$.** `QuarticRoots` calls `Solve`, so behaviour on quartics whose roots are general algebraic numbers (returned as `Root[…]` objects) depends on Mathematica's ability to simplify expressions involving them. For polynomials with rational or simple radical roots — including all examples in the paper — this works cleanly.

3. **Branch / sign issues.** The intermediate `PowerExpand` in `substUSquared` assumes the algebraic identity $\sqrt{u^2}=u$ holds (i.e. we're working at the level of formal differentials). The final closed form is correct as an antiderivative on each branch.

4. **`Integrate` as the elementary-integration oracle.** After the reduction, the integral is genuinely elementary — a rational function divided by $\sqrt{ax^2+bx+c}$. `Integrate` handles these robustly; if it ever fails to find a closed form, the issue is purely the elementary-integration step, not Goursat's reduction.
