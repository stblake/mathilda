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

---

# Correlation and template matching

## ImageCorrelate

`ImageCorrelate[image, kernel]` correlates — the kernel is **not** reflected, which is the only
difference from `ImageConvolve`. `ImageCorrelate[image, template, "NormalizedCrossCorrelation"]` is
template matching. Attributes: `Protected`.

**Correlation is convolution with the kernel reversed on both axes**, and the identity is asserted rather
than assumed:

```
ImageCorrelate[img, k] == ImageConvolve[img, Reverse[Reverse[k], 2]]
```

**bit-exactly**, because correlation is *implemented* that way — the kernel is reversed and handed to the
convolution path — so the identity holds by construction rather than by two implementations agreeing. That
also means correlation inherits **separability**: a 5×5 box is rank 1, so it costs 10 multiply-adds per
pixel rather than 25.

The first version had its own nested loop, and both consequences showed up: the identity held only to
3.6e-15 (the same products summed in a different order), and a 5×5 correlation took 3.67 ms against
scipy's 2.8 because it silently opted out of every optimisation the convolution path had accumulated.
Deriving it removed the duplicate code, the discrepancy and the deficit together.

The distinction between the two only shows on an *asymmetric* kernel, where a delta with `{{1,2,3}}` gives
`{3,2,1}` here and `{1,2,3}` convolved. Both directions are pinned, so neither can drift toward the other.

## Normalised cross-correlation

`ImageCorrelate[image, template, "NormalizedCrossCorrelation"]` is template matching. A
colour image is reduced to luminance first, since NCC compares *shape* — a property of
brightness — and combining per-channel scores would need an arbitrary rule.

### Summed-area tables

NCC needs three quantities per window: the cross term, the window sum, and the window sum
of squares. Written directly that is three sweeps of the template over every pixel. Two
identities remove the statistics from the inner loop, with `m = kw*kh`:

```
Σ (I - Ī)(T - T̄) = Σ I·T − (Σ I)·T̄
Σ (I - Ī)²       = Σ I² − (Σ I)²/m
```

so the only per-window quantities left are `Σ I` and `Σ I²`, and a summed-area table
answers each in **four lookups regardless of template size**. The cross term is then a
plain correlation, so it goes through the shared correlation path and inherits its
separable case. The tables are built through the same `clampi` the direct loop used, so
they encode exactly the edge-replicated border at no extra memory.

**One property is traded, and it is worth naming.** The direct form computed `Σ ds·dt` and
`Σ ds·ds` from the same values in the same order, so a template matched against itself gave
a peak of *exactly* 1.0. Now the numerator comes from the correlation and the variance from
the tables, so the peak is 1.0 to within a few ulp — measured at 3.3e-15. The **argmax** is
unaffected, being an integer, and that is what template matching rests on, so the suite
asserts the argmax exactly, the peak to 1e-12, and agreement with the written-out
definition to 1e-11.

A uniform window has no variance and scores 0 rather than dividing; the variance is clamped
at zero first, because cancellation can push it a hair negative and `sqrt` of that is a NaN
that would spread through `Max` and `Position` without ever looking like an error.

### Measured

512×512, construction excluded:

| Template | Mathilda | SciPy | Ratio |
|----------|---------:|------:|------:|
| 8×8 | 3.55 ms | 15.03 ms | **4.2× faster** |
| 32×32 | 3.55 ms | 146.51 ms | **41× faster** |

The reference is the identical algorithm in NumPy/SciPy (integral images plus
`signal.correlate2d`), so this compares implementations rather than methods.

Both figures are flat in template size because the cross term now goes through the
transform (below). The arc for the 32×32 case was 199.6 ms with the direct cross term,
and 21.9 ms before the summed-area tables at 8×8 — the tables removed the statistics, and
the transform then removed the cross term.

## Convolution through the transform

A direct convolution costs `w·h·kw·kh` multiply-adds, which at 512×512 with a 32×32 kernel
is 275 million and measured 197 ms. `convolve_dispatch` therefore tries three things in
order: the separable factorisation (`kw + kh` taps, unbeatable when it applies), then the
transform when a cost model says it is cheaper, then the dense loop.

**The border is the whole difficulty.** A transform gives *circular* convolution while every
filter here replicates the edge, so the wrap has to be made impossible rather than accepted.
Padding the image by replication to `PH = h + kh - 1` does exactly that: the outputs that
matter sit at indices `kh-1 … kh-2+h` of the linear result and each reads only padded rows
that exist, so no output touches a wrapped one. Rounding the transform up to a 5-smooth size
adds zeros beyond the pad, which cannot reach those indices either.

No kernel flip is needed. The direct loop computes `dst[y] = Σᵢ src[clamp(y - i + ci)]·k[i]`
with `ci = kh/2`; defining `P[a] = src[clamp(a - (kh-1) + ci)]` makes
`Σᵢ P[(y + kh - 1) - i]·k[i]` equal to it, and the left side is the linear convolution of `P`
with `k` at index `y + kh - 1`. So the answer is a shifted window of the transform's output.
The kernel's transform is computed once and reused across channels.

### The crossover was measured, not guessed

The cost model has one empirical constant, and the first value chosen for it was an order of
magnitude wrong. Forcing both paths at 512×512 with non-separable kernels:

| Kernel | Direct | Transform |
|--------|-------:|----------:|
| 3×3 | 1.05 ms | ~2.7 ms |
| 5×5 | 2.60 ms | ~2.7 ms |
| 7×7 | 5.56 ms | 2.63 ms |
| 9×9 | 10.61 ms | 2.66 ms |
| 15×15 | 35.16 ms | 2.65 ms |
| 21×21 | 78.90 ms | 2.59 ms |
| 32×32 | 198.69 ms | 2.83 ms |

The transform is **flat in the kernel** — that is the point of it — so the crossover is
wherever the dense loop passes ~2.7 ms, just under 5×5. The constant is set to switch at
7×7, the first size where the transform is a clear win rather than a wash. The initial
guessed constant switched between 15×15 and 21×21, leaving 35 ms on the table at 15×15 where
2.7 ms was available.

Consequences at 512×512: **9×9 4.0× faster, 15×15 13.3×, 21×21 30.5×, 32×32 70×**. Plain
`ImageCorrelate` with a 32×32 template went 197.3 → 2.78 ms, against SciPy's 142.4 for
`correlate2d`. Separable filters are untouched by design — `GaussianFilter` at radius 8 is
3.09 ms through the factorisation, and sending it to a transform would be slower.

Below the crossover nothing changed, which is deliberate: an identity kernel stays bit-exact,
and a transform would have cost that. Agreement with the definition written out longhand is
≤1.8e-14 across the sizes that switch, including non-square kernels in both orientations and
even extents, where a transposed pad or output crop would otherwise hide.

The branch is guarded by `USE_FFTW`; `make USE_FFTW=0` compiles it out and the dense loop
answers instead.

## ImagePad and ImageCrop

`ImagePad[image, m]` pads `m` pixels on every side; `ImagePad[image, {{left, right},
{bottom, top}}]` pads each side separately. The pair is named in Mathematica's **visual**
order — bottom before top — while the data's first row is the **top** of the image, so `top`
padding adds rows at the *start* of the array. That is the reverse of how the specification
reads, and a symmetric pad cannot detect it being wrong, so the test pads one side only and
checks which end grew.

Negative amounts crop, but may not erase the image: `ImagePad[img, -9]` on a 7×5 image
declines rather than returning something zero-sized.

`ImagePad[image, m, spec]` chooses the fill:

| `spec` | Fill |
|--------|------|
| a number | that constant value (the default is 0) |
| `"Fixed"` | replicates the edge pixel — the same boundary rule the filters use, so padding and then filtering composes with it |
| `"Reflected"` | mirrors **without** repeating the edge pixel: `{1,2,3}` padded by 1 gives `{2,1,2,3,2}`, not `{1,1,2,3,3}` |

The reflection distinction is not cosmetic. Repeating the edge doubles that sample, which
biases any subsequent average toward the border. Reflection uses a period of `2n-2`, so
padding deeper than the image itself still works — padding a 3-pixel row by 4 gives
`{1,2,3,2,1,2,3}`.

`ImageCrop[image, {w, h}]` crops to `w × h` about the centre, any odd remainder going to the
right and bottom — the same floor-division convention the kernel centres use. A crop may not
enlarge.

`ImageCrop[image]` with no size instead **trims a uniform border**, which is a different
question: how much of the frame carries no information. The border colour is read from a
corner rather than assumed black, since a scanned page's margin is white and assuming black
would trim nothing. Each edge is tested independently, so a border uniform on three sides and
not the fourth trims the three. An entirely uniform image comes back unchanged — there is no
content to keep, and a zero-sized image is not an image.

### The exact identities

Both operations are index arithmetic with no interpolation, so their composition is exact:

```
ImageCrop[ImagePad[img, m], ImageDimensions[img]] === img
```

That round trip is also the test that catches an off-by-one on either side *independently*: a
pad that adds one row too many at the top and one too few at the bottom still has the right
total size, and only the round trip notices.

A **centred** crop can only invert a **symmetric** pad. With `{{1,2},{3,4}}` the content starts
at row 4 while a centred crop begins at row 3, so the exact inverse of an asymmetric pad is a
negative pad, and that is what the suite asserts:

```
ImagePad[ImagePad[img, {{1,2},{3,4}}], {{-1,-2},{-3,-4}}] === img
```

Every fill mode is invertible the same way, since none of them touches the interior.

Both are 2-D only at present; `ImagePad` on an `Image3D` declines rather than guessing an
axis order.

### Measured

512×512 float64, one channel, construction excluded from the timing:

| Operation | Mathilda | NumPy | Ratio |
|-----------|---------:|------:|------:|
| pad 16, constant | 0.110 ms | 0.050 ms (`np.pad`) | 2.2× slower |
| pad 16, reflected | 0.221 ms | 0.059 ms (`np.pad`) | 3.7× slower |
| centred crop | 0.130 ms | 0.026 ms (slice + copy) | 5.0× slower |

The first version of the pad was 0.57 ms, and the obvious suspect — a per-pixel coordinate
map where NumPy does block copies — turned out to be worth only 0.57 → 0.45. The remaining
0.35 ms was `image_load` walking the buffer element by element through `ndt_get`, which every
filter in the subsystem pays. Making that a `memcpy` for the `"Real"` float64 case (where the
unit-scaling is the identity) took the pad to 0.110 ms and `ImageConvolve` from 1.48 to
1.14 ms as a side effect. NumPy's remaining edge is that a crop there is a *view* plus one
copy, while this returns a fresh image.

## ImagePad and ImageCrop on a volume

Both accept an `Image3D`, dispatching on the image exactly as `ImageConvolve` does.
`ImagePad[volume, m]` pads all six faces; the per-axis form takes one `{lo, hi}` pair per
axis **in the order `ImageDimensions` reports**, so

```
ImagePad[volume, {{left, right}, {bottom, top}, {first slice, last slice}}]
```

The third pair is named for what it does rather than as `{front, back}`: which end of a
volume "front" means is not something to guess, and a wrong guess here is a silent
transposition rather than an error. The **height** pair is the reversed one, matching the
2-D case — Mathematica names it `{bottom, top}` while row 1 is the top of the image, so
`top` adds rows at the start of the array. Width and depth are `{low, high}` as written.
Each convention is pinned by a test that pads one axis on one side only, since a symmetric
pad cannot tell any of them apart.

The fill modes are the planar ones and apply per axis, so `"Fixed"` in depth repeats the
first slice and `"Reflected"` in depth mirrors slices without repeating it.

`ImageCrop[volume, {w, h, d}]` crops about the centre with floor division on every axis —
the same convention as the planar crop, which is what makes the round trip exact rather
than a one-voxel shift:

```
ImageCrop[ImagePad[volume, m], ImageDimensions[volume]] === volume
```

`ImageCrop[volume]` with no size **declines**. Trimming a uniform border in three
dimensions is a genuinely different question — which faces, and a shell or a box — and
declining beats picking one silently.

### Measured

64 × 96 × 128 float64 (786,432 voxels), construction excluded:

| Operation | Mathilda | NumPy | Ratio |
|-----------|---------:|------:|------:|
| pad 8, constant | 0.378 ms | 0.217 ms (`np.pad`) | 1.7× slower |
| pad 8, reflected | 0.894 ms | 0.363 ms (`np.pad`) | 2.5× slower |
| centred crop | 0.327 ms | 0.094 ms (slice + copy) | 3.5× slower |

The first volumetric pad measured 1.42 ms — 6.5× NumPy, against 2.2× for the planar pad
whose arithmetic is the same code. The cause was `image3d_load`: the `memcpy` fast path for
a `"Real"` float64 buffer had been added to `image_load` alone, so the volumetric loader was
still walking every voxel through `ndt_go`'s dtype switch. That is this subsystem's most
repeated bug in its purest form, and the measurement is what found it — the identity tests
were all passing. Adding the same path took the pad to 0.378 ms and the crop to 0.327 ms.

## Volumetric convolution through the transform

The rank-3 path takes the same three-way dispatch: separable factorisation, then the
transform when the cost model says it wins, then the dense loop. It matters more here than
in the plane, because a `kd·kh·kw` kernel is **cubic** in the radius — 9×9×9 is 729 taps per
voxel where 9×9 is 81 per pixel.

The construction is identical to the planar one: replicate the border to
`PD × PH × PW = (d+kd-1) × (h+kh-1) × (w+kw-1)`, round each axis up to a 5-smooth length, no
kernel flip, and read the answer out of the window at `(z+kd-1, y+kh-1, x+kw-1)`. The
kernel's transform is hoisted out of the channel loop.

Memory is the thing to know: a 64 × 96 × 128 volume with a 9³ kernel transforms
72 × 108 × 144 and needs roughly **36 MB** of scratch across the two real and two complex
buffers. That is another reason the cost model has to be right rather than generous — an
unnecessary transform costs memory as well as time.

### Measured

64 × 96 × 128 float64, non-separable kernels, both paths forced:

| Kernel | Taps/voxel | Direct | Transform | SciPy `ndimage.convolve` |
|--------|-----------:|-------:|----------:|-------------------------:|
| 3³ | 27 | 11.35 ms | 11.29 ms | 7.31 ms |
| 5³ | 125 | 58.31 ms | **10.89 ms** | 48.53 ms |
| 7³ | 343 | 200.28 ms | **11.55 ms** | 180.27 ms |
| 9³ | 729 | 454.83 ms | **9.53 ms** | 421.07 ms |

The transform is flat at ~11 ms, so the crossover sits exactly at 3³ where the two paths
measure equal — and the shared cost constant fitted for rank 2 puts the switch between 3³ and
5³ with no retuning, which is what the model predicted and now what the measurement says.
Staying dense at 3³ is the right call anyway: equal speed, and none of the scratch.

Against SciPy: **4.5× faster at 5³, 15.6× at 7³, 44× at 9³**. The one remaining loss is 3³,
at 11.3 ms against 7.31 — that is the 27-tap dense loop, where SciPy's `ndimage` has a
specialised C kernel and this has a general one.

Agreement with the definition written out longhand is ≤1.8e-14 at every size, including
non-cubic extents in both orderings and even extents on every axis, which is where a
transposed pad or output window would hide — a class of mistake the volumetric paths in this
subsystem have shipped twice.

## The dense loop: interior and border

Below the transform's crossover the dense loop still runs, and it used to call the clamping
index helper on **every tap**. A pixel whose kernel reaches only rows and columns inside the
image needs no clamping at all, and for a 3×3 kernel on 512×512 that is 99.2% of pixels. So
the two cases are separated: the interior is a plain dot product over a contiguous window,
the border keeps the clamped form, and the kernel is reversed on every axis once up front so
the interior's index expressions ascend — a descending stride is not something a compiler
will vectorise.

### What it was worth, including where it was not

512×512 planar and 64 × 96 × 128 volumetric, dense sizes only (above the crossover the
transform runs instead):

| | before | after |
|---|---:|---:|
| planar 3×3 | 1.05 ms | 1.05 ms |
| planar 5×5 | 2.60 ms | 2.42 ms |
| volumetric 3³ | 11.35 ms | 7.50 ms |

**The planar split is a wash and is kept anyway.** At 3×3 the timing is unchanged across
repeated runs; at 5×5 it is about 7% better. The clamping was evidently not the bottleneck in
the plane — at nine taps per pixel the per-pixel overhead is, and the branches were
predictable enough that removing them bought nothing.

The volumetric split is worth **1.5×**, and the reason is that every tap there clamps *three*
axes: a 3×3×3 kernel does 81 clamped index computations to produce 27 multiply-adds. That
closes the last case where SciPy led — `ndimage.convolve` at 3³ measures 7.31 ms against 7.50
here, which is parity.

The planar version is kept for symmetry with the rank-3 one rather than for its timing. Two
different shapes for the same algorithm at two ranks is exactly the divergence that produced
the two volumetric bugs recorded above, and identical shapes are worth something even when
one of them is not faster.

### Why not vImage

`vImageConvolve_PlanarF` is already linked through Accelerate and is genuinely fast. It is
also **float32**. Every filter in this subsystem is float64, and the suite asserts agreement
with the written-out definition at 1e-14; single precision answers to about 1e-7. Quietly
dropping six digits of a numeric function's accuracy to make it faster is not a trade to make
invisibly in a computer algebra system. If it is ever wanted it should be an explicit opt-in
that says so in its name, not the default behind an existing head.

## Corner detection

A corner is where the gradient points in **two** independent directions, which is a statement
about the second-moment matrix of the gradient over a neighbourhood — the structure tensor:

```
M = [ ⟨Ix·Ix⟩  ⟨Ix·Iy⟩ ]        ⟨·⟩ = Gaussian-weighted average over the window
    [ ⟨Ix·Iy⟩  ⟨Iy·Iy⟩ ]
```

Reading its two eigenvalues is the whole method. Both small: flat. One large: an edge, where
the gradient has a single direction. Both large: a corner. The three distinct entries are
gradient products, each smoothed, so this composes entirely from parts already here — the
derivative stencils, and a Gaussian built as an outer product so `convolve_dispatch` factors
it and the smoothing costs `2n` taps rather than `n²`.

`CornerFilter[image]` gives the response map, `CornerFilter[image, r]` sets the window radius
(default 2), and `CornerFilter[image, r, method]` selects:

| method | response | note |
|--------|----------|------|
| `"MinimumEigenvalue"` (default) | λ_min (Shi-Tomasi) | *is* "how much does the weaker direction vary"; comparable across images |
| `"Harris"` | `det − 0.04·trace²` | cheaper (no square root); goes negative on edges, which is informative |

λ_min is computed from `√((Sxx−Syy)² + 4Sxy²)` rather than `√(trace² − 4det)`. The first is a
sum of squares and cannot go negative; the second can, through cancellation, and would give a
NaN that spreads silently.

`ImageCorners[image, r, t, d, n]` gives corner **positions**, with the window radius `r`
(default 2), the threshold `t` as a fraction of the largest response (0.05), the minimum
separation `d` in pixels (0), and the maximum number of features `n` (all). Three filters apply
**in that order**, because each removes what the others cannot:

1. **3×3 non-maximum suppression.** Alone it returns a maximum in every flat region, since a
   plateau of zeros has maxima too. A plateau yields exactly one position, by comparing strictly
   against earlier neighbours and loosely against later ones.
2. **Threshold.** Alone it returns a blob of adjacent pixels around every corner, because the
   response is smooth.
3. **Minimum separation**, which is what makes the list usable. On a noise-like 512×512 image
   the first two leave **4104** positions — every one a genuine local maximum above the cut, and
   useless as a feature set, because they arrive in clusters a pixel or two apart. Selection is
   greedy in *descending response order*: walk the sorted list and keep a position if it is at
   least `d` from everything already kept. That ordering is what makes the choice principled —
   the survivor of a cluster is its strongest member, not whichever came first in raster order.

`n` truncates last, deliberately: applied before separation it would return `n` positions out of
a single cluster.

The result is sorted strongest first, ties broken by position, so `First` is the strongest corner
and the same image always gives the same list. The greedy pass is O(kept × candidates), and **measurement says it is cheap**: at 512×512 the
raw list is 4104 corners in 6.14 ms, and `d = 10` keeps 1044 in 6.72 ms — the greedy selection
costs about 0.6 ms of that, because each candidate is rejected by the first kept neighbour it is
near and the loop exits early. A grid would make it asymptotically linear and there is no reason
to write one.

Positions carry a caveat worth reading twice if you are writing a test against them: Mathilda's
real comparison is **tolerant**, in Mathematica's way, so two responses one ulp apart are `Equal`
and neither is `Greater`. `Sort[values, Greater]` therefore cannot serve as a descending-order
oracle for a response list — it reports such elements incomparable and permutes them freely. The
suite asserts non-increasing order elementwise instead, and asserts "the top n" as *the weakest
kept is at least as strong as the strongest dropped*.

Positions are `{row, column}`, 1-based, so each indexes `ImageData` directly. This is **not**
Mathematica's `{x, y}` measured from the bottom left; the convention is stated rather than
guessed, because a silently transposed or flipped coordinate looks plausible on any test image.

### The properties that are exact

- A uniform region has no gradient, so the response is **exactly** `0.0`, not merely small.
- **A straight edge scores exactly zero** under both methods. Every gradient in the window is
  parallel, so `M` has rank 1 and both `det` and `λ_min` vanish. This is the property that
  separates a corner detector from an edge detector, and it is the one no visual inspection of
  a response map reliably shows.
- **Rotational covariance.** A quarter turn is an exact index permutation, so the response of
  the rotated image is the rotation of the response — measured agreement 1.4e-17, and four
  turns reproduce the original map under `===`. This is what catches a swapped `Ix`/`Iy`, a
  transposed smoothing pass, or a sign error in one derivative, none of which spoil the look
  of the output.
- A 24×24 checkerboard of 6-pixel blocks has interior corners at every multiple of 6, so
  `ImageCorners` returns exactly **9**.

### Measured

512×512, against scikit-image:

| | Mathilda | scikit-image | Ratio |
|---|---:|---:|---:|
| λ_min response | 5.28 ms | 10.79 ms (`corner_shi_tomasi`) | **2.0× faster** |
| Harris response | 5.33 ms | 9.08 ms (`corner_harris`) | **1.7× faster** |
| positions | 6.01 ms | 29.8 ms (`corner_shi_tomasi` + `corner_peaks`) | **5.0× faster** |

Not quite like-for-like: scikit-image is parameterised by `sigma` where this takes an integer
radius, and `corner_peaks` applies its own minimum-separation rule. The comparison is of the
same work in the same shape, not of identical outputs.

## Corner detection in a volume

`CornerFilter` and `ImageCorners` both accept an `Image3D`. In three dimensions the structure
tensor is symmetric 3×3, and **its rank says what the neighbourhood contains**:

| rank | content | eigenvalues |
|-----:|---------|-------------|
| 0 | flat | all three zero |
| 1 | a **plane** | two zero |
| 2 | an **edge** — a line where two planes meet | one zero |
| 3 | a **corner** | all three positive |

So `λ_min` is exactly zero for a plane *and* for an edge, and positive only at a true corner.
That hierarchy is the test: a detector that fires on a planar interface has not been written for
three dimensions, and no inspection of a response map shows the difference. Measured: flat gives
exactly `0.0`, a plane `< 1e-15`, an edge `9.8e-17`, and an octant corner `0.0527` at the corner
voxel.

`λ_min` comes from the closed trigonometric solution for a symmetric 3×3's eigenvalues, not an
iteration — the characteristic polynomial is a cubic with three real roots, and a per-voxel
iterative solver would need a convergence story. `r` is clamped into `[-1, 1]` before `acos`,
where rounding could otherwise produce a NaN. `"Harris"` subtracts `k·trace³`, **not** `trace²`:
`det` scales as `λ³` in three dimensions, and the two terms have to share a dimension or the
constant `k` means nothing.

Gradients are central differences taken directly on the buffer. The sign convention is
irrelevant here — the tensor holds *products* of gradients, so flipping an axis changes nothing —
which sidesteps the pre-flipped-stencil trap the planar derivative code documents. Smoothing is
three 1-D passes rather than one cubic kernel: a radius-2 Gaussian is 125 taps per voxel as a
cube and 15 as three lines, and there are six tensor entries to smooth.

`ImageCorners` on a volume returns `{slice, row, column}` triples, suppresses over the 26
neighbours, and separates in 3-D Euclidean distance. The peak-finding is **one implementation
parameterised by depth**, shared with the planar path — a plane is depth 1, where the
neighbourhood collapses to eight and the position to a pair. That is deliberate: the volumetric
paths in this subsystem have twice diverged from their planar twin by exactly one dropped detail.

A 24³ checkerboard of 6-voxel blocks returns exactly **27** corners (3³, at every multiple of 6);
a plane returns `{}`; an edge returns none.

### Measured

32 × 48 × 64 (98,304 voxels), radius 2:

| | Mathilda | scikit-image | Ratio |
|---|---:|---:|---:|
| λ_min response | 4.90 ms | 36.65 ms (`structure_tensor` + `structure_tensor_eigenvalues`) | **7.5× faster** |
| Harris response | 3.69 ms | 7.67 ms (`structure_tensor` alone) | **2.1× faster than the tensor alone** |
| positions | 4.98 ms | — | — |

Harris needs no eigenvalues, which is why it is the cheaper of the two here. Not quite
like-for-like: scikit-image is parameterised by `sigma` where this takes an integer radius.

## Named options

`CornerFilter` and `ImageCorners` take Wolfram-style trailing options as well as their positional
forms:

```
CornerFilter[image, r, Method -> "Harris"]
ImageCorners[image, MaxFeatures -> 5]
ImageCorners[image, 2, 0.05, 10., MaxFeatures -> 3]
```

`Options[CornerFilter]` reports `{Method -> "MinimumEigenvalue"}` and `Options[ImageCorners]`
reports `{MaxFeatures -> Infinity}`. Those registered defaults are the **single place** the reader
looks, so `SetOptions[ImageCorners, MaxFeatures -> 4]` takes effect with no further code in the
builtin — which is the argument for keeping defaults in the symbol registry rather than as
constants in C.

The mechanism is `options_extract` in `options.c`, and it is general: any builtin can accept a set
of named options without growing a positional tail. `ImageCorners` had reached five positional
arguments, four of them settings, which is what prompted it. Options must be **trailing**, as in
Mathematica; the scan runs from the end of the argument list so a positional argument that happens
to be a rule cannot be mistaken for an option. Duplicates resolve last-wins.

An unknown option **declines** rather than being ignored. Mathematica warns and continues; refusing
is the more conservative reading and matches this tree's rule of refusing rather than guessing — a
mistyped option name that silently does nothing is the failure it avoids.

### One thing the reader has to get right

An option entry reports whether **the call** supplied it or the value came from the defaults. That
distinction is load-bearing for any head that also accepts the setting positionally, and the first
version of this reader omitted it. The consequence was immediate and silent:
`CornerFilter[image, 2, "Harris"]` computed the *MinimumEigenvalue* response, because the default
`Method` filled the slot and the builtin could not tell that apart from an explicit option, so the
default overrode the positional argument. The default is therefore consulted only when no
positional was given.

It was caught by asserting that the option form and the positional form agree — an equivalence
worth writing down for any setting that has two spellings, since neither spelling looks wrong on
its own.

## LocalAdaptiveBinarize

A global threshold cannot binarize unevenly lit content, and that is not a tuning problem: if one
half of a page is darker than the other, **no single number** separates ink from paper in both
halves at once. `LocalAdaptiveBinarize[image, r]` compares each pixel to the mean of its own
`(2r+1)²` neighbourhood; `LocalAdaptiveBinarize[image, r, {c1, c2, c3}]` uses

```
threshold(y, x) = c1·mean + c2·stddev + c3
```

Mean alone (the default `{1, 0, 0}`) is Bradley's method; a negative `c2` is Sauvola's, tightening
the threshold where the neighbourhood is busy. Summed-area tables make the statistics **O(1) per
pixel regardless of r** — the same identity NCC uses — where a radius-16 window would otherwise be
1089 taps per pixel. With `c2 = 0` the sum-of-squares table is not built at all, since it exists
only for the deviation term. The result is typed `"Bit"`.

### The test that shows it earning its place

A checkerboard under a strong lighting ramp. No single number separates the two tones in both
halves, so a global threshold **must** collapse one half to a single value, while the local form
keeps the pattern in both. Asserted as set membership — "this half contains both values" — which is
absolute, rather than as a percentage recovered. Measured: the global result's left third is
`{0.}`; both thirds of the local result contain `{0., 1.}`.

### Mean-only is a boundary case, stated precisely

With `c2 = 0` and `c3 = 0` the threshold *is* the window mean, so a pixel can tie it — and on a
periodic synthetic image many do exactly. The summed-area mean (four lookups, then `sum/area`) and
a direct sum of the same window differ in the last bit, and a **binary** decision amplifies that
into a visible flip: 15 pixels of 224 in the suite's case, every one of them a tie. So the test
does not assert the outputs are equal; it asserts the stronger true statement that **every
disagreement is a pixel within one ulp of its own threshold**. Where the definition is numerically
determined, the implementation matches it — and with `c2` or `c1` moving the threshold off the mean,
agreement is exact.

### Measured

512×512, flat in radius as the tables promise:

| | Mathilda | scikit-image |
|---|---:|---:|
| r = 2 | 0.58 ms | 1.72 ms (`threshold_local`, mean) |
| r = 8 | 0.60 ms | 1.67 ms |
| r = 16 | 0.66 ms | 1.68 ms |
| r = 32 | 0.73 ms | 1.72 ms |

**2.6–2.9× faster** — but it was *4.2× slower* an hour before, and the reason is worth recording as
the third instance of one pattern. `bit_image_from_mask` was building 262144 `Expr` integers in 512
nested `List`s, and that dominated: global `Binarize` measured 6.56 ms where its Otsu pass and
comparison together are well under 1 ms, and the local version measured 7.2 ms — only 0.7 ms more,
despite building summed-area tables over the whole image. Emitting a packed `int64` buffer instead
took the local form to **0.58 ms** and global `Binarize` from 6.47 to **0.278 ms**, a 23× speedup of
a head this change was not even about.

The other two instances were `image_load` walking a buffer element-by-element (found when `ImagePad`
looked 10× slower than NumPy) and `image3d_load` still doing so after `image_load` was fixed (found
when the volumetric pad looked 6.5× slower). In all three the answers were correct and only the
marshalling was slow, so no test could have caught them — a benchmark did, each time.

## make check-image-packing

Three times in this subsystem an image operation ran 4× to 23× slower than its equivalent
elsewhere, with **entirely correct answers**, because the marshalling and not the algorithm was the
cost:

| where | found because |
|-------|---------------|
| `image_load` walked an NDArray element-by-element | `ImagePad` measured 0.57 ms against NumPy's 0.050 |
| `image3d_load` still walked after `image_load` was fixed | the volumetric pad measured 6.5× NumPy while the planar one had just reached 2.2× |
| `bit_image_from_mask` built 262144 `Expr` integers in nested `List`s | `LocalAdaptiveBinarize` measured 4.2× scikit-image — and the fix made global `Binarize` 23× faster |

No test in the suite could have caught any of them, because nothing was wrong with the output. A
benchmark caught each one, by accident, one at a time. So the property is now checked mechanically:
**an image-returning head hands back a packed buffer.**

Heads are read out of the image sources rather than listed in the tool, so a new one is covered the
day it is registered, and each is asked at **both ranks** — two of the three bugs were a volumetric
path failing to gain what its planar twin already had, so a planar-only gate would be half a gate.
Call shapes are discovered by trial, the same approach `nd_fastpath_sweep.py` takes and for the same
reason: a curated list of shapes only covers the heads someone remembered. It ratchets on
`KNOWN_UNPACKED`, currently empty.

It also reports, as information rather than failure, which heads accept no rank-3 shape at all —
`Binarize`, `Closing`, `Dilation`, `Erosion`, `EdgeDetect`, `MeanFilter`, `MedianFilter`, `Opening`,
`DistanceTransform`, `GradientFilter`, `DerivativeFilter`, `ColorConvert`, `ImageCorrelate`,
`ImageReflect`, `ImageRotate`, `LocalAdaptiveBinarize`. That is a genuine backlog of volumetric
coverage, listed where it can be seen instead of rediscovered.

### The gate proves it can fail

The first version of this tool reported `0 image-returning heads` and `no newly nested image heads`
and exited **0** — because a missing bracket made the generated probe a syntax error, and `run()`
returned stdout while ignoring the exit status. A gate that passes because its probe never ran is
worse than no gate: it is a green light with nothing behind it. So the exit status is checked, an
empty result is a failure rather than an empty success, and four paths are exercised directly — a
newly nested head exits 1, a head fixed since being listed exits 0 with a FIXED note, an empty probe
exits 1, and a failed probe exits 1 and keeps the generated file for inspection.

## Morphology in a volume

`Dilation`, `Erosion`, `Opening` and `Closing` accept an `Image3D` with an integer radius.

**Separability is what makes it usable.** A radius-4 box in three dimensions is 729 voxels per
output, and a max over a box is the max over lines along each axis in turn — so three van Herk
passes cost O(1) per voxel per axis and the whole operation becomes **independent of the radius**,
exactly as the planar version is. Written out directly the same filter would be cubic in `r`.

Only an integer radius is accepted at rank 3. An arbitrary 3-D structuring element is not separable
in general, so honouring it would mean the direct cubic walk — hundreds of times more work behind
the same spelling — and `Dilation[volume, element]` **declines** rather than quietly costing that.

### The properties are algebraic identities

Which is the strongest thing a morphological operator offers: they hold exactly, or the
implementation is wrong.

- **Agreement with the definition is `===`, not a tolerance.** Three van Herk passes equal the max
  over the whole cube exactly — max and min are order statistics, so unlike a sum there is no
  rounding to hide behind.
- **Idempotence** holds exactly: `Opening[Opening[v, r], r] === Opening[v, r]`. That is what makes
  an opening an opening rather than merely a smoother, and it holds only because both passes use the
  same element.
- **Duality**: eroding `f` is dilating `-f` and negating. Checked through `1 - f` so values stay in
  the unit interval, which costs one rounding — hence 1e-16 rather than exact.
- **Ordering**: `Opening[v] ≤ v ≤ Closing[v]` pointwise.
- **Monotone** in the radius, and `r = 0` is exactly the identity.
- A uniform volume is unchanged by both.

### Measured

32 × 48 × 64 (98,304 voxels):

| radius | box | Mathilda | SciPy `grey_dilation` | Ratio |
|-------:|----:|---------:|----------------------:|------:|
| 1 | 27 | 0.535 ms | 1.11 ms | **2.1× faster** |
| 2 | 125 | 0.536 ms | 1.14 ms | **2.1×** |
| 4 | 729 | 0.561 ms | 1.17 ms | **2.1×** |
| 8 | 4913 | 0.584 ms | 1.31 ms | **2.2×** |

`Opening` at r = 4 is 1.12 ms against SciPy's `grey_opening` at 2.34 ms, also 2.1×. Both
implementations are near-flat in the radius, since SciPy separates a full box too; the ratio is a
comparison of implementations, not of methods.

`make check-image-packing` now reports 37 (head, rank) pairs, all packed — the four heads added here
appear at rank 3 and the planar-only list is four shorter.

## MeanFilter and MedianFilter in a volume

Both accept an `Image3D` with an integer radius, and they are worth reading together because one is
separable and the other provably is not.

**`MeanFilter` is a box sum, done by prefix sums along each axis** — three passes, each O(1) per
voxel, so the cost is **independent of the radius**. The planar half convolves with a normalised box,
which factors into two 1-D passes but still costs `2r+1` taps per axis; a box *sum* needs no taps at
all, since differencing two prefix entries gives the window total in constant time. That is why
SciPy's `uniform_filter` is flat in `r` where a separable convolution is not — and the measurement
said so before this was rewritten: 0.42 ms at r = 1 rising to 0.94 at r = 4, against SciPy's steady
0.58.

The trade is real and named in the source. Summing only the `2r+1` window values carries error
`~k·ε`; differencing two prefix sums carries `~line·ε`, about ten times more on a 64-long line. Both
are ~1e-14 on unit-interval data, so it buys radius-independence for rounding no caller of an image
mean can observe. But it does mean **`r = 0` is no longer bit-exact by construction**, so that case
short-circuits to a copy — an identity that was silently lost when the prefix path replaced the
convolution, and that the suite caught.

**`MedianFilter` cannot do the same, and the reason is worth stating: the median is not separable.**
Taking medians along x, then y, then z gives a different filter — it has its own uses and it is not
the median of the cube. So the window really is gathered, `(2r+1)³` values per voxel, and the only
saving available is in the selection. `window_median` is therefore a **quickselect** with a
median-of-three pivot, shared by both ranks: insertion sort is the right choice at 25 elements (a
radius-2 plane) and the wrong one at 125 (a radius-2 cube), where it costs ~3900 comparisons per voxel
against quickselect's ~250. One implementation serves both ranks because it computes the same order
statistic, rather than a fast path for volumes and a slow one for planes.

### Properties

- The median matches the definition **exactly** (`===`) at both radii — an order statistic has no
  rounding to hide behind — which also pins the even-window convention, since the lower middle is
  taken rather than averaging (averaging would invent a value not in the window).
- **Every output value is an input value.** A mean invents values; a rank filter must not.
- **A lone impulse is removed completely.** In an otherwise constant volume the median returns the
  constant and *nothing else*, because the impulse is outvoted in every window containing it. A mean
  cannot do this at any radius — it spreads the impulse instead. This is the property that says why a
  median filter exists, asserted as the exact set of output values.
- `r = 0` is exactly the identity for both; a uniform volume is unchanged by both.

### Measured

32 × 48 × 64 (98,304 voxels):

| | Mathilda | SciPy | Ratio |
|---|---:|---:|---:|
| MeanFilter r = 1 | 0.550 ms | 0.58 ms (`uniform_filter`) | 1.05× faster |
| MeanFilter r = 2 | 0.558 ms | 0.57 ms | 1.02× |
| MeanFilter r = 4 | 0.557 ms | 0.60 ms | 1.08× |
| MedianFilter r = 1 | 2.81 ms | 5.28 ms (`median_filter`) | **1.9× faster** |
| MedianFilter r = 2 | 13.84 ms | 37.96 ms | **2.7× faster** |

The median grows as the window does, which is inherent — it is the one filter here with no
separable form. `make check-image-packing` now reports 39 (head, rank) pairs, all packed.

## Binarization in a volume

`Binarize` and `LocalAdaptiveBinarize` accept an `Image3D`. Both return a packed `"Bit"` volume.

`Binarize[volume]` uses Otsu over the whole volume, which needs no rank awareness at all — a
volume's histogram is a histogram, so `img_otsu` runs on the flat buffer unchanged.
`Binarize[volume, t]` is exactly a comparison against `t`.

### The 3-D window statistics reuse the box sum rather than a summed-volume table

The textbook route for `LocalAdaptiveBinarize` at rank 3 is to extend the summed-area idea one
dimension: each box sum becomes **eight lookups with alternating signs**, an inclusion–exclusion whose
sign pattern is genuinely easy to get wrong and which produces plausible-looking output when it is.

That formula is not written here. `mean3_boxsum` — three separable prefix passes, O(1) per voxel,
already checked against the definition for `MeanFilter` — gives the window mean directly, and the mean
of squares is the same call on the squared volume, from which the variance is `E[x²] − E[x]²`. Reusing
a tested routine beats hand-writing an error-prone one, and it is also cheaper in memory: two
volume-sized buffers rather than a padded table of `(D+2r+1)(H+2r+1)(W+2r+1)` doubles, which at
64×96×128 with r = 4 would be 8.4 MB per table with two tables needed.

As at rank 2, the sum-of-squares work is skipped entirely when `c2 = 0`.

### Properties

- Agreement with the longhand definition is **exact** whenever the threshold is moved off the window
  mean (`{1, -0.3, 0.02}` and `{0.9, 0, 0}` both `===`).
- Mean-only is the same boundary case as at rank 2, and confirmed to behave identically: with the
  threshold *equal* to the window mean a voxel can tie it, and the prefix-sum mean differs from a
  direct sum in the last bit. Two voxels of 1680 in the suite's case, and the assertion is that **every
  disagreement is a tie**, not that the outputs match.
- A uniform volume gives all zero; raising the offset can only turn voxels off.
- **The discriminating test at rank 3**: a checkerboard volume under a lighting ramp along x. No single
  number separates the tones at both ends, so a global threshold collapses one end to a single value
  while the local form keeps the pattern at both — asserted as set membership.

### Measured

32 × 48 × 64 (98,304 voxels):

| | Mathilda | Reference | Ratio |
|---|---:|---:|---:|
| `Binarize` (Otsu) | 0.096 ms | 0.59 ms (`skimage.threshold_otsu`) | **6.2× faster** |
| `LocalAdaptiveBinarize` r = 2 | 0.562 ms | 0.59 ms (`uniform_filter` + compare) | 1.05× |
| r = 4 | 0.529 ms | 0.61 ms | 1.15× |
| r = 8 | 0.557 ms | 0.62 ms | 1.11× |
| Sauvola form, r = 4 | 1.062 ms | 1.38 ms (two `uniform_filter`s) | 1.30× |

Flat in the radius, as the prefix-sum construction requires. The Sauvola form costs about twice the
mean-only one, which is the second box sum over the squared volume and is why it is skipped when
`c2 = 0`.
