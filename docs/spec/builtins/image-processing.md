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
