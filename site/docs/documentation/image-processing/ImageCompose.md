# ImageCompose

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageCompose[base, over] alpha-composites over onto base, centred, keeping base's size and clipping whatever falls outside. ImageCompose[base, over, {x, y}] centres the overlay at {x, y} in image coordinates -- x from the left, y from the BOTTOM. ImageCompose[base, {over, a}] scales the overlay's opacity by a. A grey image composed with a colour one produces colour: grey means the same value in every channel, so it is replicated rather than zero-padded.`**

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (7)

```mathematica
In[1]:= a = Image[Table[N[(i + j)/32], {i, 1, 16}, {j, 1, 16}], "Real"];

In[2]:= red = Image[Table[{1., 0., 0.}, {i, 1, 6}, {j, 1, 6}], "Real"];
```

```mathematica
In[3]:= ImageCompose[a, red]
Out[3]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABB0lEQVR42u3bOw6DMBQFUYMi2aWX6d2QXYRdkoLWdIjP85nSomLkm1EkplrrljqUUpLz885zzt3zOeFWCCCAABAwLh+Vck3tHD3vBpggAkAAAVBBsWtHBZkgEEAACFBBY9aOCjJBIIAAEKCCxqwdFWSCQAABIEAF7fzW9VEv4rssKsgEgQACQIAKuuDX/2nc9R7cABNEAAggACpIBcEEEQACCIAKUkEwQQSAAAKgglQQTBABIIAAqCAVBBNEAAggACcxtda2iF+gv+XcDTBBBIAAApBe/l+Q2lFBJggEEAACYleQ2lFBJggEEAACYleQ2lFBJggEEAACYleQ2lFBJggEEAAChuQPXBIcIPVlTRMAAAAASUVORK5CYII=)

```mathematica
In[4]:= ImageDimensions[ImageCompose[a, red]]
Out[4]= {16, 16}

In[5]:= ImageChannels[ImageCompose[a, red]]
Out[5]= 3
```

```mathematica
In[6]:= ImageCompose[a, red, {4, 4}]
Out[6]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABAElEQVR42u3bOwqAMBQF0acISZmd6y50l7oB7cRP3pkyWAiDlynM0Frb44Raazi/77yUcno+Bl6FAAIIAAF5mVTKM7Vz9bwvwAQRAAIIgArqu3ZUkAkCAQSAABWUs3ZUkAkCAQSAABWUs3ZUkAkCAQSAABWUs3ZUkAkCAQSAABWUs3ZUkAkCAQSAABWUs3ZUkAkCAQSAABWUs3ZUkAkCAQSAgL9U0Lptn3rRZZ5VEEwQASCAADxQQV+j1/+XfAEmiAAQQABUkAqCCSIABBAAFaSCYIIIAAEEQAWpIJggAkAAAVBBKggmiAAQQADuqqCrG+i93slSQSCAABBAACIi4gCNwRuaUhXx7gAAAABJRU5ErkJggg==)

```mathematica
In[7]:= ImageCompose[a, {red, 0.4}]
Out[7]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABJ0lEQVR42u3csQ6DIBRAUdpgZOSLXP3/ratf0g5dYVNQOHd8MR284XnTRF85528okFIK5ufN13Utzt8BXSGAAAJAwLxEldKmdmrXOwFWEAEggACooLFrRwVZQSCAABCgguasHRVkBYEAAkCACpqzdlSQFQQCCAABKujPfhzF+bIsl85jjMX5Z9tUkBUEAggAASqowdO/V+3Uru91H5wAK4gAEEAAVFCb2lFBIIAAEEAAnlVBZ9WOCgIBBIAAAnDPCrq6dmq/r4KsIBBAAAhQQW2e/r1qx39BIIAAEEAA2lRQ7Q30Xu9k3W3uBFhBBIAAAvD0Cur1vR0VBCuIABBAANpUkNpRQVYQCCAABIxdQWpHBVlBIIAAEDB2BakdFWQFgQACQMCU/ACOxBzt1Z2GygAAAABJRU5ErkJggg==)

### Applications (3)

```mathematica
In[8]:= a = Image[Table[N[(i + j)/32], {i, 1, 32}, {j, 1, 32}], "Real"];
```

An edge map laid over the image it came from

```mathematica
In[9]:= ImageCompose[a, {EdgeDetect[a], 0.6}]
Out[9]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAC8klEQVR42u2dvaoiMRiG4yoWXoJ3YO0NDAiCYCdY2QmCiqBYiIWICIqVjQiCgrUgWlgJVlbegPdilS3cDVl2RjM/Z+YcfN7qK/LHDO9D8mWSiSWTSSn+KJFI/A1FPB5X8Xq9VnGj0bAt46euXqbZbIY2Bqd+Tfpyijebzdu+dP0SKFLxAiJWLJVKSTcWM7F8mO1EVddtO6vVSsWtVgsHgCD0dA/YCQc7/zx0rS4OAEEgKBC76QsNt+24XaA59eXUjr44MinvNAaT8Tthp9PpqHg+n+MAEIReI8gPdnRLup11mOSFosKOSV0n7Ohler0eDgBB6JkLOhwO0o3lo8KOSYo4Kuwsl8u32HGKcQAI+nAEpdNpCXbCxc5kMsEBIAh5S0eDHW/Y0ePRaIQDQBAKZ0fs07Azm81UPBgMbPsdDoc4AAShYHfEmO283nzXy0ynUxwAgpD/HTGw8xpBTmX6/T4OAEFICCFELJPJSLDjDzvj8dj24eqzHXbEQBByWCuAnbCxA4JAEFKzoNvtJsM6Rgp2nmq32zgABKEngrLZrPSajgY7/ncScQAIIh391jLf4UzWT8FOt9tV8WKxsC1Tr9dxAAhC/yPIJC8Edrxhx+lZ4QAQ9OELsfv9LoO4ngvseBsDDgBBH44gy7Ik2AkXO9wXBIKQq6+jwY45dmq1moq32y3paBCEgv06Guz4xw4IAkFILcQKhYIM4lJTsOPtGeIAEMTX0WAnQOxUKhUV73Y723ZKpRIOAEHI22UdYOc1doTBMdXT6YQDQBAyRxDYMcfOfr8XfB0NgtCX3RcEdvxjBwSBIKRc4nTOy22OCOyYP7d8Po8DQBASQggRq1QqEuwEh51isaji8/nMLAgEodcIejwebw9omPyez+QuZRP9dOzoulwutnUty8IBIAg9EVStVmXQ2OFn0BzQAEHIQzoa7IQ/fhwAgj58FiSllG7+ga7/jPgrbBvmV8rlclnFx+PRc24nl8up+Hq9siMGgpCxfgPmQkvA89pd+AAAAABJRU5ErkJggg==)

A blurred copy blended halfway: the classic soft-focus composite

```mathematica
In[10]:= ImageCompose[a, {GaussianFilter[a, 3], 0.5}]
Out[10]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABQklEQVR42u3dsQ6CMBhF4WIM8AB9//dkqosDA0SI8EPLd6YOLnhzTgxD7XLOJX0ZhiEtncdxbPIc+Yx936clXgmXYoCLec81nGuylqO1z9d43vKM/5zXssMACcJigvae7/yr5g5nBkgQfiZoi8K1ZOQOGWSABKG6BNX+3okBEoRDEjR/dyE752WHARKESxIkOwyQIAQnSHYYIEEITpDsMECCEJwg2WGABCE4QbLDAAlCcIJkhwEShOAEyQ4DJAgS1Ex2GGAANJkgBsAAEiQ7DDAA2k4QA2CA5hMkOwyQIDw8QQyAASRIdhhgALSdIAbAAC5uDb6yDAwwAIJvTQQDDIDgBIEBBkBwgsAAAyA4QWCAARCcIDDAAAhOEBhgAAQnCAwwAM6hm6appB3/gQ4GGAAHJqiUUnwNDDAADGAAGOBxfADUrDaW3eaixgAAAABJRU5ErkJggg==)

### Properties & Relations (5)

```mathematica
In[11]:= a = Image[Table[N[(i + j)/32], {i, 1, 16}, {j, 1, 16}], "Real"];

In[12]:= red = Image[Table[{1., 0., 0.}, {i, 1, 6}, {j, 1, 6}], "Real"];
```

The size is the base's, whichever way round the two are given

```mathematica
In[13]:= {ImageDimensions[ImageCompose[a, red]], ImageDimensions[ImageCompose[red, a]]}
Out[13]= {{16, 16}, {6, 6}}
```

Outside the overlay, the grey base is replicated across all three channels

```mathematica
In[14]:= Module[{d = ImageData[ImageCompose[a, red]]}, d[[1, 1, 1]] === d[[1, 1, 2]] && d[[1, 1, 2]] === d[[1, 1, 3]]]
Out[14]= True
```

At zero opacity the overlay contributes nothing, even where it covers

```mathematica
In[15]:= Module[{d = ImageData[ImageCompose[a, {red, 0.}]]}, d[[8, 8, 1]] === d[[8, 8, 2]]]
Out[15]= True
```

## Implementation notes

- `Protected`.
- The result keeps the **base's** size and clips whatever falls outside: composition is "draw on
  this", not "make something bigger".
- A grey image composed with a colour one produces colour. Grey means the same value in every
  channel, so it is **replicated**, never zero-padded — padding would turn a grey pixel red.
- The result carries alpha only if the base did: compositing onto an opaque image gives an opaque
  image.

**Attributes:** `Protected`.

## References

- Source: [`src/imagecompose.c`](https://github.com/stblake/mathilda/blob/main/src/imagecompose.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
