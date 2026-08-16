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

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAADBklEQVR42u2dT+spURjHj0sWXoWXYGunlLJTVspCKeXfQpTYSMSGWEgRZaeUsvAClIWlnRfiDcxduPc0C2POMePnuj7f1aPmOcPo++k5z5k54/H7/Yb4I5/P9zcUXq9XxovFQsaFQuHuMU5y3RrnFbm68XK5tD3XbDaT8S+B3ir+gDfLEwgEjGft5pZtnYzzrlzdcczYwQEgCEn3gJ2fx06lUsEBIAipI0hlgmOVa56M6NrfPL45V2VCZ5Wrcl6r8c2TLKvjVbAzGo1wAAhCjxGkix23eiyvwI5KrkpvR+X3WmHHfEy9XscBIAjdnGq2bbFYBDtPYGc6nd69uCrnwgEg6MsRZIUd3UkK2LmpVquBIBCEXteOBjvOsdPpdHAACEJCCCE8wWDQADuvx47VdcMBIIgqCOy4iJ1eryfjdrvNRAwEIfV29L+2gP6J1Y6gHQ2CkA6CwI67vR0QBIKQFoLAjnPsdLtdrdxqtYoDQBASQgjhCYVCBthxDzuNRkPGw+HQtreGA0AQ7Wih044GO4+xYzURY0UMBKG7VdD5fDbe3Y7+n7AzHo+1riEOAEFfjqBwOGyAnZ/FTi6XwwEgCN0QdLlcDCZZ7mGnVCrJeD6fUwWBIKS+KA92nFc7KovytKNBEHppO/rbsJPP54VOOzqdTuMAEISe2y8I7DzGzmq1ss3dbDY4AAQh9f2CwI5z7NALAkHIlf2CwI5z7CSTSRwAgpAQQghPIpEwwI572EmlUjLebre0o0EQcnfLMrDzGDuCu6NBEHIbQWBHHTu73c42Nx6P4wAQhG4Tsev1avuAhsrbSMEO7WgQhJ5AUDqdfvoBDbDz3HWLRqM4AAShx/sFmSsfFeluamoltx6OsNK7sGOOj8cjDgBB6FYFZbNZw+5lxOVy+eNepvwpL4PGASCIWxNVSiWw86LvjwNA0JcjyPxBpdoZDAYybjabtnZrtVoy7vf7QmdT08lkYtvbsVqty2QyMl6v11q3C+73e9vjY7GYjA+Hg4wjkYiMT6cT7WgQhB7qNzrWFfNstz+1AAAAAElFTkSuQmCC)

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
