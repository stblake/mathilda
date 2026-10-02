# Histogram

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Histogram[data, opts...]`**

Bins the numeric values in data and draws a frequency histogram. Bin count defaults to Sturges' rule: ceil(Log2\[n\]) + 1.

**`Histogram[data, k, opts...]`**

k equal-width bins.

**`Histogram[data, {step}, opts...]`**

Bins of width step.

**`Histogram[data, {min, max, step}, opts...]`**

Explicit range and width. Options: ChartStyle   color/style list cycling through bins BarSpacing   gap fraction of bin width (default 0.2) Standard Graphics options pass through.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Histogram[Table[RandomReal[], {200}]]
Out[1]= -Graphics-

In[2]:= data = Table[RandomReal[], {200}]
Out[2]= {0.140681, 0.719171, 0.778469, 0.0741672, 0.813947, 0.580171, 0.013126, 0.754838, 0.466119, 0.69305, 0.169174, 0.693866, 0.692065, 0.703363, 0.792567, 0.786629, 0.865939, 0.0184981, 0.224943, 0.947115, 0.825416, 0.872362, 0.425489, 0.711869, 0.882917, 0.2937, 0.957468, 0.814742, 0.952457, 0.736825, 0.507261, 0.108332, 0.600101, 0.72911, 0.921463, 0.563866, 0.30962, 0.84094, 0.0390951, 0.596077, 0.71277, 0.981742, 0.242259, 0.405348, 0.432131, 0.146109, 0.842136, 0.771742, 0.329693, 0.0130025, 0.159663, 0.0559964, 0.467622, 0.749642, 0.0197631, 0.604011, 0.206575, 0.0286502, 0.346594, 0.930716, 0.345432, 0.0202961, 0.461868, 0.217948, 0.0192531, 0.239615, 0.549122, 0.43497, 0.825891, 0.734764, 0.0863096, 0.906252, 0.596247, 0.192763, 0.239217, 0.72976, 0.881986, 0.704707, 0.815473, 0.541946, 0.496483, 0.673814, 0.856206, 0.841128, 0.753937, 0.241638, 0.682385, 0.523032, 0.320656, 0.728542, 0.924473, 0.913413, 0.0807655, 0.642171, 0.726602, 0.802352, 0.875943, 0.975244, 0.393971, 0.324169, 0.379947, 0.381772, 0.18552, 0.779192, 0.0800715, 0.846092, 0.865819, 0.568104, 0.31528, 0.0117362, 0.0222622, 0.412189, 0.582518, 0.729162, 0.859468, 0.956868, 0.691664, 0.928095, 0.450906, 0.976306, 0.949257, 0.768445, 0.762611, 0.146943, 0.560973, 0.0514557, 0.406394, 0.667696, 0.135431, 0.439167, 0.986205, 0.616763, 0.636469, 0.63966, 0.248824, 0.729922, 0.0165414, 0.439217, 0.940262, 0.106709, 0.290777, 0.299902, 0.189051, 0.914778, 0.493506, 0.742124, 0.342366, 0.878464, 0.49417, 0.577868, 0.751581, 0.200443, 0.774339, 0.372416, 0.549911, 0.147261, 0.992123, 0.817696, 0.748853, 0.231703, 0.388255, 0.528721, 0.599479, 0.160558, 0.747748, 0.648368, 0.0302633, 0.499011, 0.563099, 0.769536, 0.799137, 0.11899, 0.184662, 0.970667, 0.724645, 0.715689, 0.785371, 0.586788, 0.0478026, 0.214353, 0.600657, 0.724945, 0.625272, 0.0148207, 0.598829, 0.913496, 0.618059, 0.740913, 0.529887, 0.00678126, 0.845692, 0.595929, 0.759282, 0.306366, 0.9679, 0.635907, 0.573828, 0.248266, 0.688067, 0.914912}

In[3]:= Histogram[data, 20]
Out[3]= -Graphics-

In[4]:= Histogram[data, {0.1}]
Out[4]= -Graphics-

In[5]:= Histogram[data, {0, 1, 0.05}]
Out[5]= -Graphics-
```

## Algorithm

barchart.c — BarChart[data, opts...] and Histogram[data, opts...]

BarChart renders a vertical bar chart from explicit heights. Histogram bins numeric data and renders a frequency histogram. Both return Graphics[...] objects auto-displayed by the REPL.

BarChart[{v1,...,vn}, opts...]

```text
  n bars at x = 1..n with heights v1..vn, coloured via ChartStyle or
  the default palette.
```

BarChart[{{v1,...}, {w1,...}, ...}, opts...]

```text
  Multiple grouped datasets, each dataset in a distinct palette colour.
```

Histogram[data, opts...]

```text
Histogram[data, k, opts...]          k equal-width bins
Histogram[data, {step}, opts...]     bins of given width
```

Histogram[data, {min,max,step}, opts...] explicit range + width

Options (both):

```text
  ChartStyle   — color/style list cycling through bars
  ChartLabels  — label list for x-axis ticks
  BarSpacing   — gap fraction of bar width (default 0.2)
  PlotLabel, standard Graphics options pass through 
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [ChartStyle](../../other-advanced/ChartStyle/), [BarSpacing](../../other-advanced/BarSpacing/)

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
