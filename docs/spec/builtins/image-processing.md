# Image Processing

The raster image subsystem. This page currently covers the **representation and its
accessors**; filtering (`ImageConvolve`, `GaussianFilter`, `Binarize`) and `Image3D` build on
these and land separately.

## Image

`Image[data]` is a raster image, normalising to the canonical `Image[data, type]`.
Attributes: `Protected`.

- `Image[data]` — type inferred from the values
- `Image[data, "Bit" | "Byte" | "Real"]` — type stated, and validated against the data

**The canonical form is `Image[data, type]`, which is also real Wolfram syntax.** Normalising to
it is what makes validity *decidable*: a builtin that returns `NULL` leaves its expression alone,
so a valid `Image[data]` and a nonsensical `Image["hello"]` would otherwise be
indistinguishable — both merely unevaluated. A canonical two-argument form is exactly what passed
validation, and `ImageQ` tests for it.

**Data layout, and the transposition that is the whole trap.** `data` is a rectangular
`height × width` array, indexed `data[[y, x]]` with rows running *down* the image; a colour image
is `height × width × channels`, interleaved. `ImageDimensions` reports `{width, height}` —
**transposed relative to this**. That is Mathematica's convention and the single most common
source of silently-wrong image code, so the tests use a non-square image throughout: a square one
cannot tell the two apart.

**Type inference looks at the values and nothing else**, so it is predictable from the data
alone:

| data | inferred type |
|---|---|
| all integers, values within `{0, 1}` | `"Bit"` |
| all integers, values within `0..255` | `"Byte"` |
| anything else | `"Real"` |

So `Image[{{0, 1}, {1, 0}}]` is a bit image while `Image[{{0., 1.}, {1., 0.}}]` is a real one —
the same numbers written differently, and a distinction a caller can rely on.

A **stated** type is refused where the data does not fit it: `Image[{{0, 300}}, "Byte"]` stays
unevaluated rather than reinterpreting 300, which would corrupt every later scaling.

**Ragged data declines** rather than being padded or truncated. A ragged array is not an image,
and every filter downstream indexes it as rectangular — so a clear refusal here beats an
out-of-bounds read somewhere far away.

## ImageQ

`ImageQ[expr]` gives `True` for a valid canonical image, `False` otherwise. Attributes:
`Protected`.

This is *how validity is tested*, because `Head` cannot do it: malformed input to `Image` stays
unevaluated, so both a valid image and a refused one have head `Image`.

## ImageDimensions

`ImageDimensions[image]` gives `{width, height}`. Attributes: `Protected`.

Transposed relative to `ImageData`, which is `height × width`.

## ImageChannels

`ImageChannels[image]` gives the number of colour channels — `1` for grey, otherwise the length
of each pixel's value list (`3` for RGB, `4` with alpha). Attributes: `Protected`.

A pixel list of a different length in one place is not a channel count, so that declines rather
than reporting the first pixel's length as if it were the image's.

## ImageType

`ImageType[image]` gives `"Bit"`, `"Byte"` or `"Real"`. Attributes: `Protected`.

The type is not decoration: it fixes the **range** of a stored value, which is what makes
`ImageData`'s scaling well defined.

## ImageData

`ImageData[image]` gives the pixels as reals in `[0, 1]`, scaling out the image's type.
`ImageData[image, type]` gives the stored values unscaled. Attributes: `Protected`.

A `"Byte"` 255 comes back as **exactly** `1.0` and a 0 as exactly `0.0`; the tests assert these
as equalities rather than tolerances, so a wrong divisor (254, or 256) fails outright.

```
In[1]:= ImageData[Image[{{0, 128, 255}}]]
Out[1]= {{0.0, 0.501961, 1.0}}

In[2]:= ImageData[Image[{{0, 128, 255}}], "Byte"]
Out[2]= {{0, 128, 255}}
```

Storing the original values and scaling **on read** keeps `Image[…]` round-trippable — what went
in is what `FullForm` shows — and makes the scaling an assertion about one function rather than a
property spread across two. Real data passes through untouched, *including* values outside
`[0, 1]`: storing faithfully beats clamping silently, which would destroy data the caller may
want back.

`ImageData[image, type]` accepts only the image's **own** type. Converting between types has its
own rounding decisions, and an accessor should not make them silently.

## Performance

Accessors are **O(height)**, not O(pixels). The first implementation routed every accessor
through the full validator, which made `ImageDimensions` cost **0.59 ms** on a 512×512 image —
asking how wide an image is touched all 262,144 pixels. Since the canonical form can only be
produced by `Image[]`, which validates every pixel once, the accessors re-check only
**rectangularity** (O(height) row-length comparisons). That is ~1000× faster:

| | before | after |
|---|---|---|
| `ImageDimensions` (512×512) | 0.59 ms | 0.58 µs |
| `ImageQ` (512×512) | 0.62 ms | 0.64 µs |

Rectangularity is still checked, and a test pins it: a hand-typed
`Image[{{1,2},{3}}, "Bit"]` must not pass as an image, because every filter downstream will index
it as rectangular.

`ImageData` on a 512×512 byte image is ~4.8 ms, dominated by building 262,144 `Expr` reals — the
cost of materialising an array as expressions, and the reason filters will operate on the stored
data rather than round-tripping through `ImageData`.

---

# Filtering

## ImageConvolve

`ImageConvolve[image, kernel]` convolves `image` with a rank-2 numeric kernel. Attributes:
`Protected`.

**This is true convolution: the kernel is reflected before summing.** That distinction is not
pedantry, and it is unusually easy to get wrong, because on any *symmetric* kernel — a Gaussian, a
box — convolution and correlation agree exactly. Every smoothing filter anyone tries first looks
right either way. It shows up only on an asymmetric kernel, where the two answers are mirror
images:

```
In[1]:= ImageData[ImageConvolve[Image[{{0., 1., 0.}}], {{1, 2, 3}}]]
Out[1]= {{1.0, 2.0, 3.0}}          (* correlation would give {{3., 2., 1.}} *)
```

Mathematica draws the same line — `ImageConvolve` reflects, `ImageCorrelate` does not — and a test
pins it with exactly that case.

**Padding is `"Fixed"`**: out-of-range reads clamp to the nearest edge pixel, replicating the
border. This is Mathematica's default, and it is the right one for smoothing — zero padding would
darken every edge, which looks exactly like a real vignetting bug. Clamping also makes an exact
property available: a constant image convolved with a kernel summing to 1 comes back as the *same*
constant everywhere, border included. A test asserts that, and it is what would fail if the
padding were ever changed.

The result is always a `"Real"` image of the same dimensions; each colour channel is convolved
independently. A Gaussian of bytes is not a byte, and rounding back into the input type would
discard precision the caller never asked to lose.

## GaussianMatrix

`GaussianMatrix[r]` gives a `(2r+1) × (2r+1)` Gaussian normalised to sum 1.
`GaussianMatrix[{r, sigma}]` states the standard deviation. Attributes: `Protected`.

`sigma` defaults to `r/2`, Mathematica's convention: it puts the kernel's edge at two standard
deviations, where the Gaussian has fallen to about 13% of its peak, so truncating there loses
little.

**Normalisation divides by the realised sum, not the analytic `2 pi sigma^2`.** The analytic
constant is correct only for an infinite kernel; using it on a truncated one leaves the sum
slightly under 1, which darkens an image a little on every pass — invisible once, obvious after
fifty.

## BoxMatrix

`BoxMatrix[r]` gives a `(2r+1) × (2r+1)` matrix of 1s. Attributes: `Protected`.

**Not normalised**, matching Mathematica — so `ImageConvolve[img, BoxMatrix[1]]` is nine times too
bright, and the normalised version is a mean filter. Kept faithful rather than helpfully rescaled,
since a caller using `BoxMatrix` in arithmetic needs the ones.

## GaussianFilter

`GaussianFilter[image, r]` blurs `image` with a Gaussian of radius `r`. Attributes: `Protected`.

Exactly `ImageConvolve[image, GaussianMatrix[r]]` — Mathematica documents the two as equal, so it
is implemented by building the same matrix and calling the same convolution, and a test asserts the
identity. Two independent implementations of one identity is how the identity quietly stops
holding.

## Filtering performance, measured

Verified against `scipy.ndimage.convolve` (the right baseline: it reflects the kernel like
`ImageConvolve`, and `mode='nearest'` replicates the edge exactly like `"Fixed"`). The two agree to
**six significant figures** on sampled pixels and on the image total — 2050.0623 either side.

512×512 grey, normalised Gaussian:

| radius | taps | Mathilda | scipy | ratio |
|---|---|---|---|---|
| r=0 | 1 | 4.6 ms | — | marshalling floor |
| r=1 | 9 | 5.3 ms | — | |
| r=2 | 25 | **6.3 ms** | 3.2 ms | **1.97×** |
| r=4 | 81 | 14.4 ms | — | |

**The `Expr` round-trip dominates at small radii, not the inner loop.** A single-tap convolution
still costs 4.6 ms, which is purely loading 262,144 pixels into a buffer and rebuilding the result
as expressions; the convolution arithmetic at r=2 is only ~1.7 ms of the 6.3. Arithmetic scales
linearly in taps (1.7 ms at 25 taps, ~9.8 ms at 81), which is what the direct form should do.

That measurement redirects the optimisation target. The inner loop is not the problem; the
expression marshalling is. In priority order: keep images in packed buffers so a filter chain
does not round-trip through `Expr` at every step; exploit **separability** for Gaussians
(`kw + kh` taps instead of `kw × kh`, a 2.5× win at r=2 and 4.5× at r=4); and only then reach for
Accelerate's `vImage`, which is already linked. Nothing here is GPU-accelerated, and at ~1.7 ms of
arithmetic per megapixel-scale image the transfer cost to a GPU would exceed the saving —
`vImage` is the honest next step, not Metal.

The benchmark pair lives in `benchmarks/63-image-convolve/`, including a deliberate
`r=0` row that measures the marshalling floor on its own.

---

# Thresholding and colour

## FindThreshold

`FindThreshold[image]` gives a threshold separating the image into two classes, by **Otsu's
method**. Attributes: `Protected`.

Otsu maximises the **between-class** variance

```
sigma_b^2(t) = w0(t) w1(t) (mu0(t) - mu1(t))^2
```

which is algebraically the same as minimising the weighted *within*-class variance — so one cheap
pass optimises the thing you actually want (how well separated the two groups are) without ever
computing a within-class variance. Both weights and both means update **incrementally** as `t`
advances one bin, so the search is O(bins) after the histogram; recomputing the means per
candidate is the obvious implementation and is quadratic.

A colour image is reduced to luminance first. Values are binned over `[0, 1]`, and anything
outside — which a `"Real"` image may legitimately hold, since `Image` stores faithfully rather than
clamping — is clamped *into* the histogram rather than dropped, so an out-of-range pixel still
votes for the extreme it belongs to instead of vanishing from the statistics.

**Returns unevaluated for an image whose pixels are all identical.** That is one cluster; no
threshold splits it into two, and inventing one would be a fiction.

**On comparing against other libraries.** The returned value is the winning bin's **upper edge**,
binned over `[0, 1]`. scikit-image returns the bin **centre**, binned over `[min, max]` of the
data. Both are legitimate, and they differ by a bin or two while agreeing on *which split is
best* — so the test verifies the defining property instead of matching a number: between-class
variance is recomputed from the pixel values in Mathilda itself, sharing no code with the
implementation, and the returned threshold must maximise it. On a 48×48 ramp it does, exactly.

## Binarize

`Binarize[image]` thresholds by Otsu; `Binarize[image, t]` thresholds at `t`. Gives a `"Bit"`
image. Attributes: `Protected`.

**A pixel strictly above the threshold becomes 1**, so a pixel exactly at it becomes 0. That
matters and is pinned by a test: "above" and "at or above" differ on exactly the pixels a threshold
was chosen to sit between.

```
In[1]:= ImageData[Binarize[Image[{{0.4, 0.5, 0.6}}], 0.5], "Bit"]
Out[1]= {{0, 0, 1}}
```

On a bimodal image there is a ground truth, and Otsu recovers it exactly — every pixel on the right
side, not merely most:

```
In[2]:= ImageData[Binarize[Image[{{0.2,0.2,0.8,0.8},{0.2,0.2,0.8,0.8}}]], "Bit"]
Out[2]= {{0, 0, 1, 1}, {0, 0, 1, 1}}
```

The result is built as a `"Bit"` image directly rather than as `"Real"`: it is 0/1 by construction,
and typing it `"Real"` would mean a caller could no longer tell it was binary.

## ColorConvert

`ColorConvert[image, "Grayscale"]` converts by **Rec. 601 luminance**,
`0.299 R + 0.587 G + 0.114 B`. Attributes: `Protected`.

Weighted rather than averaged because the eye is far more sensitive to green than to blue. An
unweighted mean would give 1/3 for each of pure red, green and blue, putting a saturated blue and a
saturated green at the same brightness — and for thresholding specifically it would put pure red
and pure blue on the same side of any threshold, when perceptually they are far apart. The test
asserts the three weights as exact values for exactly this reason.

**Only `"Grayscale"` is accepted.** The other spaces Mathematica supports (LAB, HSB, XYZ, …) each
carry their own white point and transfer-function decisions; accepting the name while doing
something approximate would be worse than declining it.

---

# Geometry

## ImageResize

`ImageResize[image, {w, h}]` resizes to `w × h`; `ImageResize[image, w]` gives width `w` with the
height following to preserve the aspect ratio. `Resampling -> "Nearest" | "Bilinear" | "Average"`
selects the method. Attributes: `Protected`.

**Aliasing is the whole problem with downsampling, and it is invisible until it is catastrophic.**
Picking every other pixel out of a photograph looks fine, so nearest-neighbour reduction passes
every casual test. Hand it a fine checkerboard and it returns a flat field — the pattern is sampled
at exactly the frequency that annihilates it:

```
In[1]:= chk = Image[{{0.,1.,0.,1.},{1.,0.,1.,0.},{0.,1.,0.,1.},{1.,0.,1.,0.}}];

In[2]:= ImageData[ImageResize[chk, {2, 2}]]                      (* area averaging *)
Out[2]= {{0.5, 0.5}, {0.5, 0.5}}

In[3]:= ImageData[ImageResize[chk, {2, 2}, Resampling -> "Nearest"]]
Out[3]= {{0.0, 0.0}, {0.0, 0.0}}                                  (* pattern gone *)
```

Nyquist requires every frequency above half the *new* sampling rate to be removed **before**
resampling, and no interpolation afterwards can restore what point-sampling discarded. So
`Automatic` uses **area averaging whenever either axis shrinks** — a box prefilter and a resample in
one pass. It is not the best possible antialiasing filter (a windowed sinc has better stopband
behaviour) but it is exactly right for integer reduction factors, cheap, and impossible to get
subtly wrong. Enlarging has no frequencies to remove, so bilinear is used there; area averaging on
an enlargement would degenerate to nearest, since each destination pixel would fall inside a single
source pixel.

**Fractional coverage, not integer blocks.** Each source pixel is weighted by how much of it the
destination pixel actually overlaps, so a 3 → 2 reduction is as correct as 4 → 2. Restricting the
area path to integer factors would have been simpler and would have quietly fallen back to
something worse on the sizes people actually ask for.

**Centre-aligned coordinates.** A destination pixel `i` covers source `[i·s, (i+1)·s)` with its
centre at `(i+0.5)·s`, so the bilinear map is `sx = (i+0.5)·s − 0.5`. The naive `sx = i·s` is the
classic half-pixel shift: correct at 1:1, drifting the image half a pixel at every other scale, and
asymmetric — the left edge gains a border the right does not. A test *only* at 1:1 would miss it
entirely, which is why the checkerboard and fractional rows exist alongside the identity ones.

An unknown resampling name **declines** rather than falling back to a default, so a typo cannot look
like it worked. Sizes must be positive integers: rounding `10.5` silently would make the call mean
something the caller did not say.

### Measured

Verified against skimage (`downscale_local_mean` for the area path — exactly a block mean at an
integer factor, which is what area averaging reduces to; `resize(order=1, anti_aliasing=False)` for
bilinear). The checkerboard reduction agrees exactly, `[[0.5, 0.5], [0.5, 0.5]]` either side.

| operation | Mathilda | skimage | ratio |
|---|---|---|---|
| 512→256, area | **1.71 ms** | 2.0 ms | **0.86× (faster)** |
| 512→256, nearest | 1.67 ms | — | |
| 512→1024, bilinear | 17.8 ms | 14.3 ms | 1.24× |

**Area averaging preserves the mean exactly** — 0.497985 before and after, to every digit — because
every source pixel contributes with equal total weight. That is a conservation law, so the benchmark
carries it as a *check* rather than a timing, and it is what fails if the coverage weights are ever
normalised wrongly.

The two ratios together confirm the earlier diagnosis rather than adding a new one. Mathilda is
*faster* on the reduction, where the output is a quarter of the input and there are few `Expr`s to
build; it is slower on the enlargement, where 17.8 ms for 1M output pixels matches the ~4.6 ms per
262k marshalling rate measured for convolution almost exactly. **The resampling arithmetic is not
the cost — building the output expressions is**, and that cost scales with the *output* pixel count.

Benchmark pair in `benchmarks/64-image-resize/`.

---

# Storage: machine buffers

Computed images store their pixels as a **machine buffer** (a rank-2 or rank-3 `NDArray`) rather
than nested expressions. This is invisible through the API — `ImageData` always answers with a
`List`, whatever the storage — and it is what three separate benchmarks had been pointing at.

| operation, 512×512 | nested `Expr` | buffer | speedup | baseline |
|---|---|---|---|---|
| marshalling floor (r=0) | 4.62 ms | **0.98 ms** | 4.7× | — |
| `ImageConvolve` r=2 | 7.18 ms | **2.95 ms** | 2.4× | scipy 3.2 ms |
| `ImageConvolve` r=4 | 14.6 ms | **11.2 ms** | 1.3× | — |
| `ImageResize` 512→1024 | 16.4 ms | **2.76 ms** | 5.9× | skimage 14.3 ms |
| `ImageResize` 512→256 | 1.72 ms | **0.63 ms** | 2.7× | skimage 2.0 ms |
| chain of three `GaussianFilter` | 17.3 ms | **4.43 ms** | 3.9× | — |

Mathilda is now **faster than scipy** on convolution (2.95 ms against 3.2) and **5.2× faster than
skimage** on bilinear enlargement. Pixel values are unchanged — still agreeing with
`scipy.ndimage.convolve` to six significant figures.

**The visible-`NDArray` surface is required, and this is the subtle part.** Mathilda has two array
surfaces: a *packed List*, which looks like an ordinary `List`, and a *visible* `NDArray`. The
obvious choice is the packed List, since an image's pixels are conceptually a list — and it does not
work. `eval.c` has a **post-gate**: when a node comes to rest still holding a packed List, it
materialises it into expressions **unconditionally**, with no `packed_aware` check. That is right for
an ordinary head, where a resting buffer means some fast path declined it and an inert `Mod[buffer]`
would behave differently from `Mod[{1., 2., 3.}]`. It is fatal for a *container*: `Image[…]` coming
to rest holding its pixels is the entire point. A packed-List image was materialised on **every
evaluation**, which is why adding the image heads to `AWARE` changed nothing measurable.

The gate never touches a visible `NDArray`, so that is the surface image storage uses — both for
filter output and for pixels handed in as a packed `List` (which `Table` produces).

Storage is buffer-backed **uniformly, at every size**; there is no threshold at which the
representation changes. That matters for testing: it means no class of image that only large-input
tests would reach. The data a caller passes to `Image[data]` directly is kept as given.

---

# Separable kernels

`ImageConvolve` detects a **separable** (rank-1) kernel and runs two one-dimensional passes instead
of one two-dimensional one, turning `kw × kh` multiply-adds per pixel into `kw + kh` — at radius 8,
34 instead of 289.

| radius | taps 2-D → 1-D | Mathilda | scipy `convolve` (2-D) | scipy `gaussian_filter` (separable) |
|---|---|---|---|---|
| r=2 | 25 → 10 | **1.48 ms** | 2.5 ms | 1.6 ms |
| r=4 | 81 → 18 | **1.97 ms** | 10.5 ms | 1.9 ms |
| r=8 | 289 → 34 | **3.44 ms** | 53.2 ms | 3.3 ms |

So Mathilda now **matches scipy's dedicated separable routine** while being 5–15× faster than its
general 2-D convolve. The difference is not that scipy is slow: `scipy.ndimage.convolve` simply does
not detect separability, so you have to know to reach for `gaussian_filter` instead. Here **any**
rank-1 kernel gets it automatically — a Gaussian, a box, a smoothing-times-derivative outer product.

**Separability lives in `ImageConvolve`, not only in `GaussianFilter`, and that is deliberate.** Put
in one of the two, it would have made the documented
`GaussianFilter[i, r] == ImageConvolve[i, GaussianMatrix[r]]` identity *approximate*, since two 1-D
passes sum in a different order from one 2-D pass. Both routing through the same path keeps the
identity **bit-exact**, and a test asserts it.

The two-pass result equals the direct form for a real reason rather than by luck: `"Fixed"` padding is
itself separable — the direct read `src[clamp(y)][clamp(x)]` clamps the axes independently, which is
exactly what doing one axis then the other does. Only summation order differs.

**The tolerance is the whole risk, and it is tight.** Treating a non-separable kernel as separable
does not make the answer slightly wrong; it computes a completely different filter. Detection is a
relative check at `1e-12` against the kernel's own magnitude: a Gaussian factorises to ~1e-16
relative, so real separable kernels pass with room to spare, while anything genuinely rank 2 fails
long before it could be mistaken. A test convolves with `{{1,0},{0,1}}` — rank 2 — and pins the
hand-computed direct answer `{{1,0,0},{0,1,0},{0,0,0}}`, which no rank-1 approximation produces.

The factorisation pivots on the **largest-magnitude** entry rather than on `K[[1,1]]`, so a kernel
with a zero corner — which any derivative kernel has — factors instead of dividing by zero. A test
covers `{{0,0,0},{1,2,1},{0,0,0}}` for that.

Verification is independent rather than self-referential: an outer product `{{1,2},{2,4}}` must give
the same result as convolving with `{{1,2}}` and then `{{1},{2}}` — two genuine 1-D convolutions
reached by a different route through the code.

---

# Derivatives and gradients

## DerivativeFilter

`DerivativeFilter[image, {n, m}]` gives the `n`-th derivative down the rows and the `m`-th across
the columns, each order 0–2. Attributes: `Protected`.

The kernel is a separable outer product of 1-D stencils:

| order | stencil | on what it is exact |
|---|---|---|
| 0 | `{1, 2, 1}/4` | smoothing; preserves a constant *and* a linear ramp exactly |
| 1 | central difference | `f(x) = c x` → exactly `c` |
| 2 | `{1, -2, 1}` | `f(x) = c x²` → exactly `2c` |

So `{0, 1}` is Sobel-x and `{1, 0}` is Sobel-y. The stencils are **normalised**, unlike the raw
integer Sobel kernels, which report a gradient eight times the true slope — harmless when only the
*ranking* of edges matters, and wrong for anything that reads the number.

**The sign convention, and the bug it hid.** `ImageConvolve` reflects its kernel, so a central
difference written in natural reading order `{-½, 0, ½}` computes the **negated** derivative. A
derivative filter is really a *correlation*, which is why every CV library defines Sobel with
correlate semantics. The stencil is therefore pre-flipped to `{½, 0, -½}`, and convolution then
yields the true derivative — positive where brightness increases to the right, matching
`scipy.ndimage.correlate` to the digit.

That was a live bug, and it is worth recording *why the obvious tests missed it*: the gradient
**magnitude** squares the sign away, so it looked perfectly correct, and so would any edge-detection
result, since only the ranking matters there. Only an exact *signed* value could see it — "the
derivative of a ramp of slope ⅛ is exactly +⅛". An absolute-property test on the magnitude would have
shipped it.

Both stencils are separable, and the first has a **zero in the middle** — precisely the case the
separable path's largest-magnitude pivot exists for. The kernels are built as full 2-D matrices and
handed to the same dispatcher every other filter uses, so the factorisation is *re-derived and
verified* rather than trusted because the author knew it was separable.

## GradientFilter

`GradientFilter[image]` gives the gradient magnitude `Sqrt[dx² + dy²]`. Attributes: `Protected`.

**The magnitude rather than `|dx| + |dy|`, because it is rotation invariant.** An edge at 45° reports
the same strength as one at 0°; the absolute sum would report it `Sqrt[2]` times stronger and so bias
every downstream threshold by orientation. A test asserts exactly that: a diagonal ramp `x + y` of
slope ⅛ has magnitude `0.125 Sqrt[2]`.

Colour is reduced to luminance **first** and differentiated once, not differentiated per channel and
combined — per-channel gradients would need a combining rule (max? sum? norm?) and every choice is
arbitrary, where the gradient of brightness needs none.

### Measured

| operation, 512×512 | Mathilda | scipy |
|---|---|---|
| `DerivativeFilter[{0,1}]` | 1.58 ms | — |
| `GradientFilter` | **2.25 ms** | 2.3 ms |

Parity, and the values agree exactly — gradient total 26191.5 and maximum 0.522601 on both sides.

---

# Edge detection

## EdgeDetect

`EdgeDetect[image]` finds edges by the **Canny** algorithm, giving a `"Bit"` image.
`EdgeDetect[image, r]` sets the Gaussian smoothing radius (default 2; `0` means none).
`EdgeDetect[image, r, t]` sets the high threshold. Attributes: `Protected`.

Four stages, each fixing a specific failure of the one before:

1. **Smooth** — a derivative amplifies noise (differencing doubles noise amplitude while a real edge
   keeps its step), so gradients are taken of a blurred image.
2. **Gradient** — magnitude and direction from the normalised Sobel pair.
3. **Non-maximum suppression** — a magnitude ridge is several pixels wide, so thresholding alone
   gives a thick band. NMS keeps only pixels maximal *along the gradient direction*.
4. **Hysteresis** — one threshold either breaks long edges where they weaken or admits noise
   everywhere. Two thresholds plus 8-connectivity keeps a weak pixel only if reachable from a strong
   one.

**Thinning to exactly one pixel depends on an asymmetric comparison**, and this is the subtle part. A
clean step does *not* give a single-pixel gradient peak: with the central difference, a step at column
`k` responds `0.5` at **both** `k−1` and `k` — an even plateau. Testing `mag >= both neighbours` keeps
both and gives a two-pixel edge. So the test is `mag > backward && mag >= forward`: ties resolve to the
lower-index side, deterministically, and the ridge thins to one.

```
In[1]:= Map[Total, ImageData[EdgeDetect[step, 0], "Bit"]]
Out[1]= {1, 1, 1, 1, 1, 1, 1, 1}
```

The high threshold defaults to **Otsu on the *suppressed* magnitude**, not the raw one: after
thinning, the histogram really is edge-against-non-edge, which is the two-class problem Otsu solves.
On the raw magnitude it would be dominated by the wide ridge flanks. Low is `0.4 × high`, the
conventional ratio — stated because it is a choice, not a derivation.

**NMS is insensitive to the gradient sign**, worth noting after the `DerivativeFilter` sign bug: it
compares both neighbours along the direction, and direction and direction+180° select the same pair.
So that bug would not have been caught here either — more reason the exact signed value had to be
asserted where it was.

### Two deliberate differences from scikit-image

Measured against `skimage.feature.canny` on the same 8×8 step, and they disagree in two ways that are
worth stating rather than smoothing over:

| | Mathilda | skimage |
|---|---|---|
| edge width on a perfect step | **1** | 2 |
| border rows | kept | discarded |
| 512×512, r≈1–2 | **9.1 ms** | 19.9 ms |

The width difference is the plateau tie: skimage compares symmetrically, so on an *exactly* even
ridge it keeps both pixels. On real photographs the plateau is rare — values differ slightly — so both
give one-pixel edges in practice; the synthetic step is precisely the degenerate case that exposes the
tie rule. Mathilda's is the stricter reading of "one pixel wide".

For the border, off-the-edge counts as **zero magnitude**, so a border pixel is kept if it beats the
one neighbour it has. Clamping instead would compare a pixel against itself and keep every border
pixel unconditionally; discarding the border, as skimage does, loses real edges running along it.

Mathilda is ~2.2× faster. Edge *fractions* are not comparable between the two, since the thresholds are
chosen by different rules, and reporting them side by side would imply a comparison that is not there.

---

# Volumes

## Image3D

`Image3D[data]` is a volumetric image, normalising to `Image3D[data, type]`. Attributes:
`Protected`.

Data is **depth × height × width** — slices outermost, indexed `data[[z, y, x]]` — or
`depth × height × width × channels` for colour. `ImageDimensions` reports
**{width, height, depth}**: *fully reversed*.

That reversal is the 2-D transposition trap made worse. With three axes there are **six** possible
orderings, and a cubic test volume validates none of them — so every test here uses distinct
extents (2 slices of 3 rows of 4 columns), and asserts *both* directions, since either alone would
pass with two axes swapped:

```
In[1]:= v = Image3D[Table[N[(x + 10 y + 100 z)/1000.], {z, 2}, {y, 3}, {x, 4}]];

In[2]:= {ImageDimensions[v], Dimensions[ImageData[v]]}
Out[2]= {{4, 3, 2}, {2, 3, 4}}

In[3]:= Part[ImageData[v], 2, 3, 4]        (* x=4, y=3, z=2 *)
Out[3]= 0.234
```

**A volume is not a plane.** `ImageQ` is `False` for an `Image3D` and `Image3DQ` is `False` for an
`Image` — deliberately, because every filter written for a plane would otherwise accept a volume
silently and index it wrongly.

Everything else is shared with `Image`: type inference (`"Bit"`/`"Byte"`/`"Real"` from the values), the
refusal of a stated type the data does not fit, ragged rejection, the canonical form as a fixed point,
and `ImageData`'s scaling — a byte volume returns 255 as exactly `1.0`, the same code path as a byte
plane. `ImageDimensions`, `ImageChannels`, `ImageType` and `ImageData` all accept either rank.

A **colour volume is rank 4**, which is what forced `ImageData`'s nested rebuild to become a general
recursion over dims rather than the unrolled two-and-a-bit levels the plane case used. Unrolling a
fourth level would have been the point at which the pattern should have been a recursion from the
start.

Volumetric *operations* — 3-D convolution and filtering — are the next step and are not here yet.
Separability pays even better in three dimensions: `kw + kh + kd` taps instead of `kw · kh · kd`, so
27 becomes 9 at radius 1 and 729 becomes 27 at radius 4.

## Volumetric convolution

`ImageConvolve[volume, kernel]` takes a rank-3 kernel; `GaussianFilter[volume, r]` builds a 3-D
Gaussian. Dispatch is on the **image**, not the kernel: a rank-3 kernel handed to a plane is a
mistake worth declining, not something to reinterpret as a stack of 2-D kernels.

**Separability matters more in three dimensions than in two, and by a widening margin.** A rank-1
kernel costs `kw + kh + kd` taps instead of `kw · kh · kd`: radius 1 goes from 27 to 9, radius 4 from
729 to 27. In 2-D, skipping separability costs a factor of the radius; here it costs its **square**.

64³ volume (262,144 voxels):

| radius | taps 3-D → 1-D | Mathilda | scipy `convolve` (3-D) | scipy `gaussian_filter` |
|---|---|---|---|---|
| r=1 | 27 → 9 | **2.28 ms** | 3.0 ms | 2.4 ms |
| r=2 | 125 → 15 | **2.11 ms** | 18.0 ms | 2.5 ms |
| r=4 | 729 → 27 | **2.92 ms** | 155.1 ms | 2.9 ms |

At radius 4 that is **53× faster** than scipy's general 3-D convolve, and level with its dedicated
separable routine. Again the gap is not scipy being slow — `ndimage.convolve` does not detect
separability, so a caller must know to reach for `gaussian_filter`. Any rank-1 kernel gets it here
automatically. The mean is preserved exactly (0.498134 either side of an r=1 filter).

The rank-3 factorisation follows the rank-2 one: pivot on the largest-magnitude entry so a
zero-cornered kernel still factors, take the three axis-lines through the pivot as candidate factors,
then **verify every entry** against their product at a tight relative tolerance. A test checks it
independently — a separable 3×3×3 kernel must give exactly what three successive 1-D convolutions
give — and a non-separable rank-3 kernel is confirmed to run on the direct path.

The z axis gets its own reflection test, since the 2-D rows cannot reach it: a delta in a 3-slice
column with a kernel varying only in z gives `{1, 2, 3}` under convolution and would give `{3, 2, 1}`
under correlation.

---

# Morphology

`Dilation[image, r]`, `Erosion[image, r]`, `Opening[image, r]`, `Closing[image, r]`, each also taking
an explicit structuring element in place of the radius. Attributes: `Protected`.

**Flat morphology**: only the element's **support** — its nonzero positions — enters the maximum or
minimum, never its values. That is what keeps `Dilation[img, BoxMatrix[1]]`, `Dilation[img, 1]` and
`Dilation[img, {{5,5,5},{5,5,5},{5,5,5}}]` the same operation, and all three are asserted equal.

## The algebraic laws, and why they are the right tests

Morphology has exact laws, and each fails for a *different* bug:

| law | fails for |
|---|---|
| `Erosion[f,k] == 1 - Dilation[1-f,k]` | a swapped min/max, or a padding rule that is not self-dual |
| `Erosion <= Opening <= f <= Closing <= Dilation` pointwise | a mis-centred element |
| `Opening[Opening[f]] == Opening[f]` | a mis-composed pair |

All three hold **exactly** (max error 0), and the ordering chain is one assertion covering four
operators *and* their relationship to the original — something no per-operator check would catch, since
a shifted element leaves each operator individually plausible.

**Padding replicates the border**, the same rule the convolutions use, and that is what makes the laws
hold *at the edges*. Zero padding would let a dilation at the border see black that is not there,
breaking `Dilation >= f` on the boundary; and it is not self-dual, so duality would fail there too.

The cleanest statement of what dilation *is*: a single bright pixel spreads to exactly the element's
footprint, so the output **is** the element. A test pins that pixel pattern.
## Measured: constant in the radius

The 1-D max and min use **van Herk–Gil-Werman**: three comparisons per pixel *whatever the element's
width*, so the operation no longer depends on the radius at all.

A window of width `k` cannot span three blocks of width `k`. So cut the line into blocks of exactly `k`,
build a prefix maximum forwards and a suffix maximum backwards within each, and any window straddles
exactly two adjacent blocks — its maximum is `max(suffix at its start, prefix at its end)`. Two lookups
and one comparison, plus one pass each to build the arrays.

512×512:

| radius | Mathilda | scipy `grey_dilation` | before van Herk |
|---|---|---|---|
| r=1 | 2.38 ms | 2.6 ms | 1.51 ms |
| r=2 | **1.98 ms** | 2.4 ms | 2.04 ms |
| r=4 | **2.01 ms** | 2.4 ms | 2.04 ms |
| r=8 | **2.02 ms** | 2.4 ms | 3.75 ms |
| r=16 | **2.03 ms** | 2.3 ms | 6.52 ms |
| r=32 | **2.03 ms** | 2.3 ms | — |
| `Opening` r=2 | **3.73 ms** | 4.8 ms | 2.74 ms |

Flat from r=2 to r=32, and at or ahead of scipy at every radius — the crossover that used to sit at
r≈4–6 is gone, and r=16 improved 3.2×. `Opening` at r=2 regressed slightly (2.74 → 3.73 ms) from the
extra scratch passes: that is the trade, a constant cost everywhere instead of a cheap small case and an
expensive large one.

**It is exact, not approximate**, because max and min are associative and idempotent — splitting a window
at a block boundary and recombining loses nothing. That has a testing consequence: the fast path must
agree **bit-exactly** with a direct computation, so the test compares it against a reference written in
Mathilda sharing no code with the C.

Comparing against `Dilation` with an all-ones matrix would **not** be a check — an all-ones element *is* a
full rectangle, so both sides take van Herk and agree with themselves. The reference instead uses `Span`
with `Max`/`Min` over the clamped index range, which is exactly what replicate padding means. A radius
larger than the image is tested separately, that being where the last short block makes a window span the
whole padded line.

## MorphologicalComponents

`MorphologicalComponents[image]` labels the connected components of the foreground.
`[image, t]` takes pixels above `t` as foreground (default 0). `CornerNeighbors -> False` selects
4-connectivity instead of the default 8. Attributes: `Protected`.

**Returns an integer matrix, not an `Image`** — a deliberate deviation from Mathematica. Here an
`Image` would be actively wrong: type inference would call a label array of 1..12 a `"Byte"` image, and
`ImageData` would then divide every label by 255. Labels are indices, not brightnesses, and the one
thing that must not happen to them is being scaled.

**Connectivity is the discriminating property**, and the only one that can tell the two rules apart:

```
In[1]:= MorphologicalComponents[Image[{{1., 0.}, {0., 1.}}]]
Out[1]= {{1, 0}, {0, 1}}                                        (* 8-conn: one component *)

In[2]:= MorphologicalComponents[Image[{{1., 0.}, {0., 1.}}], CornerNeighbors -> False]
Out[2]= {{1, 0}, {0, 2}}                                        (* 4-conn: two *)
```

Every other property — background staying 0, labels contiguous, a single blob labelling 1 — holds under
*either* rule, so a wrong default would pass all of them.

**Two passes, doing genuinely different jobs.** The first walks in raster order and can only see
already-visited neighbours (W, NW, N, NE), so a **U-shaped** region gets *different* labels on its two
arms — nothing has connected them yet — and only the base reveals they are one component. Union-find
records the equivalence; the second pass applies it. A one-pass implementation returns two components
there and looks perfectly reasonable on every convex shape, so a U is tested explicitly under
4-connectivity (where a diagonal cannot join the arms instead).

The second pass also **relabels to 1..k in raster order of first appearance**. Without it the labels
would be whatever the first pass allocated, with gaps where two provisional labels were later merged —
so `Max` would not be the component count and no label pattern would be assertable.

Dilation can only *merge* components, never split them, so the count cannot increase: four isolated
corners become one blob at radius 1, which is an absolute relationship between two features rather than
a property of either alone.

### Measured, and the cost is not the algorithm

512×512 binarised, against `scipy.ndimage.label`: both find **41 components**, exactly. Timing is
5.20 ms against scipy's 1.2 ms — 4.3× slower.

The arithmetic says where it goes. The `Expr` marshalling rate measured for convolution was ~4.6 ms per
262,144 elements, and this builds exactly 262,144 integer `Expr`s for its output matrix. So the
labelling itself is ~0.6 ms — **faster than scipy's 1.2** — and the output construction is the other
4.6. That is the same finding as the buffer-storage change, now confirmed a fourth time, and the same
remedy: return a packed integer array instead of nested expressions. It is not done here because the
result is a `List` rather than an `Image`, so it has no container to rest inside and would face the
post-gate directly — worth doing, and worth doing carefully rather than at the end of an iteration.

---

# Rank filters

`MedianFilter[image, r]` and `MeanFilter[image, r]` over a `(2r+1)²` neighbourhood. Attributes:
`Protected`.

## The median is the one operator here that is not separable

Worth stating plainly after a week of leaning on separability. A **sum**, a **maximum** and a
**minimum** all decompose — the sum over a rectangle is the sum over rows of the sums over columns, and
likewise for max and min — because all three are associative, commutative reductions that **ignore how
values are grouped**. A median does not: it depends on a value's **rank within the whole window**, and
grouping destroys rank information.

```
{{1, 2, 9}, {3, 4, 5}, {6, 7, 8}}    true median of all nine  = 5
                                      median of row medians {2,4,7} = 4
```

The separable version is fast, plausible, and simply wrong. A test pins `5` on exactly that window, so
an implementation that took the shortcut fails rather than passing quietly.

`MeanFilter`, by contrast, **is** a convolution with a normalised box, so it is implemented as one and a
test asserts the identity against `ImageConvolve` — two implementations of one identity is how the
identity quietly stops holding. Being a full rectangle it gets the separable path.

## Why a median at all

A median removes an isolated outlier **exactly**; a Gaussian or mean only attenuates and smears it. One
bright pixel in a constant field leaves *nothing* behind under the median — the whole result is the
background value — where the mean leaves a visible bump. Both halves are asserted, because "it smooths"
is true of the mean too and proves nothing. On salt-and-pepper noise that is the difference between
clean and merely blurred.

For an even window the **lower middle** is taken rather than averaging the two, so the output is always
one of the inputs — averaging would invent a value not present in the window, which for a rank filter is
precisely what a caller does not want. A test checks every output value is drawn from the input.

## Measured

512×512, against `scipy.ndimage`:

| | Mathilda | scipy |
|---|---|---|
| `MedianFilter` r=1 (9 values) | **3.13 ms** | 6.8 ms |
| `MedianFilter` r=2 (25 values) | **10.27 ms** | 15.3 ms |
| `MeanFilter` r=2 | **1.49 ms** | 1.6 ms |

Faster than scipy at both median radii. Insertion sort on the window is the reason: for the radii people
actually use — 1 to 3, so 9 to 49 values — a small contiguous array beats a histogram or a running
median on constant factors, and it is exact and obviously correct on doubles. Both implementations grow
with the window, so at large radii a histogram median (8-bit data) or a running median would win; that
is the same shape of limit as van Herk for morphology, and the same honest position — a known ceiling
with a known remedy, not numbers that are good everywhere.

---

# DistanceTransform

`DistanceTransform[image]` replaces each pixel by its **exact** Euclidean distance to the nearest
background pixel; background pixels are 0, so the value rises toward a blob's interior.
`[image, t]` takes pixels above `t` as foreground. Attributes: `Protected`.

**Exact, not the classic chamfer approximation**, and that choice is what makes the tests equalities. A
two-pass chamfer transform propagates integer step costs and cannot represent `Sqrt[2]`, so diagonal
distances come out a few percent wrong — invisible on a picture, and it would have forced a tolerance
where an equality is available.

The discriminating test is a **3-4-5 triangle**: one background pixel, and the pixel three across and
four down must read exactly `5`. Chamfer gives about 5.03. A 5-12-13 triangle is checked too, so 3-4-5
cannot be passing by coincidence, along with `Sqrt[2]` at the diagonal neighbour — the value chamfer
fundamentally cannot produce.

Uses Felzenszwalb and Huttenlocher's method: per axis it computes
`D(x) = min over y of ((x−y)² + f(y))`, the **lower envelope of parabolas**. Every parabola has the same
curvature, so any two intersect exactly once and the envelope is built in a single sweep maintaining a
stack of still-visible parabolas — O(n) per row, no sorting.

**Separability is exact here, unlike the median's**, and the reason is worth naming: squared Euclidean
distance is a **sum** over the axes, so minimising it decomposes per axis. A median does not decompose
because rank is not a sum. The same word covers an exact factorisation in one case and a wrong shortcut
in the other, and which it is depends on whether the reduced quantity is additive. The square root is
taken once at the end — per pass it would be wrong, not merely slower.

### Measured

512×512 on an identical explicit mask, against `scipy.ndimage.distance_transform_edt` (also exact):

| | Mathilda | scipy |
|---|---|---|
| time | **2.86 ms** | 7.3 ms |
| max distance | 6 | 6 |
| total | 323895 | 323895 |

2.6× faster with values agreeing exactly, which is the check that matters — both compute the same exact
transform, so any disagreement would be a bug rather than a convention difference.

A further test ties it to morphology: surviving `k` erosions is equivalent to distance `>= k+1`. The
`+1` is the substance — a 3×3 erosion removes every pixel with a background neighbour, so one erosion
already means distance 2. Both `k=1` and `k=2` are asserted so the relation is pinned rather than a
coincidence at one value.

## Volumetric resampling

`ImageResize[volume, {w, h, d}]` and `ImageResize[volume, w]` (the other two extents following the
aspect ratio), with `Resampling -> "Nearest" | "Bilinear" | "Average"`. Dispatch is on the **image**, so a
two-element spec handed to a volume declines rather than being guessed at as a plane.

The 2-D reasoning carries over unchanged, and so does the discriminator. A 2×2×2-periodic pattern halved
on every axis is sampled at exactly the frequency that annihilates it:

```
In[1]:= chk = Image3D[Table[N[Mod[x + y + z, 2]], {z, 4}, {y, 4}, {x, 4}]];

In[2]:= Union[Flatten[ImageData[ImageResize[chk, {2, 2, 2}]]]]
Out[2]= {0.5}                                              (* area averaging: the mean *)

In[3]:= Union[Flatten[ImageData[ImageResize[chk, {2, 2, 2}, Resampling -> "Nearest"]]]]
Out[3]= {0.0}                                              (* pattern gone *)
```

Fractional coverage in three axes — a source voxel's weight is the product of its three per-axis
overlaps — so an integer factor gives an exact block mean over eight voxels and a non-integer one is
handled correctly rather than approximated. Area averaging **preserves the mean** exactly at an integer
factor, a conservation law in 3-D as in 2-D.

**The axis order is the harder part, and it is not the arithmetic.** The spec is `{width, height, depth}`
while the storage is `depth × height × width` — fully reversed — so a resize must reverse the spec before
indexing. Both a non-cubic **source** and a non-cubic **target** are needed to test it: with three axes
there are six ways to get it wrong and a cube hides every one. The tests resize a 4×3×2 volume to
`{2, 6, 4}` and assert *both* the reported dimensions and the data shape (`{4, 6, 2}`), since either alone
would pass with axes swapped.

### Measured

64³ volume:

| operation | Mathilda | reference |
|---|---|---|
| 64³ → 32³, area | **0.83 ms** | skimage `downscale_local_mean` 0.8 ms |
| 64³ → 128³, trilinear | **8.08 ms** | scipy `zoom(order=1)` 36.7 ms |

Parity on the reduction and **4.5× faster** on the enlargement.

---

# Levels and tone adjustment

## ImageLevels

`ImageLevels[image]` gives `{{level, count}, …}` — the histogram as **data, not a plot**. Mathematica's
`ImageHistogram` draws a graphic and `ImageLevels` gives the counts; in a computer algebra system the
counts are the useful half, so that is what this provides and `Histogram` over the result draws the
picture. `ImageLevels[image, n]` uses `n` bins. Accepts volumes. Attributes: `Protected`.

**The counts sum to the pixel count, exactly.** That is the property worth asserting: every pixel lands
in exactly one bin, so a total that disagrees means a bin boundary is wrong or a pixel was dropped.

A `"Bit"` image uses its **2** natural levels and `"Byte"` its 256, because those *are* its distinct
values — binning them into anything else would invent structure. Only `"Real"` has no natural set, and is
binned into 256 over `[0, 1]`. Levels are reported on `ImageData`'s unit scale, so a level can be compared
against a pixel value without rescaling.

```
In[1]:= ImageLevels[Image[{{0, 1}, {1, 0}}]]
Out[1]= {{0.0, 2}, {1.0, 2}}
```

## ImageAdjust

`ImageAdjust[image]` stretches to the full range. `ImageAdjust[image, {c, b}]` and `[image, {c, b, g}]`
apply contrast, brightness and gamma. Accepts volumes. Attributes: `Protected`.

The default stretch puts the darkest pixel on **exactly 0** and the brightest on **exactly 1**, and two
exact properties follow:

- **Idempotence** — a second stretch is the identity. This is the row that catches an off-by-one in the
  range, because a slightly wrong divisor still produces a plausible-looking contrast curve; only a
  second pass reveals it.
- **Monotonicity** — a stretch may not reorder pixels.

**A constant image has no range to stretch and comes back unchanged.** Dividing by zero is not the
answer, and mapping the single value to either end would be arbitrary.

**The parametric curve is Mathilda's documented choice, not a claim of bit-compatibility** with
Mathematica, whose exact formula is not published in a form worth guessing at. It is stated so a caller
can reproduce it:

```
v' = (v − 1/2)(1 + c) + 1/2 + b      clipped to [0, 1],  then  v'^(1/g)
```

Contrast **pivots about mid-grey**, so a contrast change does not also shift brightness — mid-grey is
fixed by any contrast value, and a test pins that (`ImageAdjust[Image[{{0.5}}], {3., 0.}]` is `0.5`). The
clip precedes gamma because a negative base has no real power, and a test covers a value driven past 1 by
brightness. `{0, 0, 1}` is the exact identity.

### Measured

512×512:

| operation | Mathilda | reference |
|---|---|---|
| `ImageLevels` (256 bins) | **0.66 ms** | numpy `histogram` 1.6 ms |
| `ImageAdjust` stretch | **0.76 ms** | skimage `rescale_intensity` 7.7 ms |
| `ImageAdjust {c,b,g}` | 1.49 ms | skimage `adjust_gamma` 1.3 ms |

The stretch is 10× faster than skimage's `rescale_intensity`, which does more work per call (dtype
negotiation and range inference); the gamma path is at parity, as two `pow` loops over the same buffer
should be.
