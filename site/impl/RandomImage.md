---
source: src/imageio.c
---
**Algorithm.** `builtin_random_image` fills a buffer with uniform noise. The max value
(default 1), size (`read_size`: a single `n` means `{n, n}`, default `150 × 150`) and a
`ColorSpace -> "RGB"` option (three independent channels, vs. the one-channel `"Grayscale"`
default) are parsed through `options_extract`. Each sample is `random_uniform_01() * max` —
drawn from the **same** stream as `RandomReal`, so `SeedRandom` makes a random image
reproducible; a private generator would have quietly made this the one random builtin that
ignores the seed. The samples are scaled by `max` but **not** clamped: the caller asked for
that range, a `"Real"` image may legitimately hold values above 1, and clamping is `Export`'s
job.

`RandomImage` exists because a filter is most honestly judged on a noise field — a smoothing
radius means nothing on a checkerboard and everything on noise — and before it the only way to
get one was a deterministic `Mod` expression masquerading as random.

**Data structures.** A fresh `w · h · channels` buffer, wrapped by `image_build_real` as a
packed `"Real"` image (every image head returns a packed buffer — `make check-image-packing`).

**Complexity / limits.** `O(w · h · channels)` draws from the global RNG. Dimensions must be
positive integers (a symbolic size stays unevaluated rather than becoming a guess); an
unsupported colour space declines rather than defaulting to grey.
