# DerivativeFilter

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DerivativeFilter[image, {n, m}] gives the n-th derivative down the rows and the m-th across the columns, each order from 0 to 2. The kernel is a separable outer product of 1-D stencils: order 0 is the smoothing {1,2,1}/4, order 1 the central difference {-1,0,1}/2, order 2 the second difference {1,-2,1}. So {0,1} is Sobel-x and {1,0} is Sobel-y. The stencils are NORMALISED, unlike the raw integer Sobel kernels, which report a gradient eight times the true slope -- harmless when only the ranking of edges matters, and wrong for anything that reads the number. On f(x) = c x the first derivative gives exactly c. The result is a "Real" image.`**

## Examples (27)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (6)

```mathematica
In[1]:= ImageData[DerivativeFilter[Image[{{0., 0., 1., 1.}, {0., 0., 1., 1.}}], {0, 1}]]
Out[1]= {{0.0, 0.5, 0.5, 0.0}, {0.0, 0.5, 0.5, 0.0}}

In[2]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[3]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];
```

```mathematica
In[4]:= DerivativeFilter[ramp, {0, 1}]
Out[4]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAqklEQVR42u3RoREAIAwEQUAl/RcbCQ1EZgazJ1/+7sy8qykiluaqqnY/rvkbAAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAAMAFAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAA0EwPNBsE124pZU4AAAAASUVORK5CYII=)

```mathematica
In[5]:= DerivativeFilter[ramp, {1, 0}]
Out[5]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAm0lEQVR42u3RAQ0AAAjDsIN/zyCDhHQS1koy0VltAQAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIwIcWi54BvzXkqDQAAAAASUVORK5CYII=)

```mathematica
In[6]:= ImageDimensions[DerivativeFilter[chk, {0, 1}]]
Out[6]= {16, 16}
```

### Scope (13)

```mathematica
In[7]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[8]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[9]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[10]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[11]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];
```

```mathematica
In[12]:= DerivativeFilter[chk, {0, 1}]
Out[12]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAA50lEQVR42u3awQ2CQBBA0Y+xsC1tS6MzbUAOmCAY3r9K8PDCZAOzVK8+NOdsT67/7vpHOjUAAAAIwH17bv2wruuhf+z+ngAjSAAAANCpLW28C9rbGOPQU8S/398TYAQJAAABACAAAAQAgAB04S9i3u385v6eACMIgAAAEAAAAgBAAAAoe0HZC5IRBEAAAAgAAAEAIAAAlL2g7AXJCAIgAAAEAIAAABAAAMpeUPaCBACAAAAQAAACAEAAACh7QdkLkhEEQAAACAAAAQAgAACUvaDsBRlBAgBAAPJFzCnIKcgIEgAAAnCf3tqWa8FOlOOvAAAAAElFTkSuQmCC)

```mathematica
In[13]:= DerivativeFilter[chk, {1, 1}]
Out[13]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAzUlEQVR42u3awQnAIBAAwTWkcDtPmhBFmH1LHhm8h9yovnSsxy8AAEAAAOhM76oPzTmddwOMIAEAIAAABACAAAAQgO58C/K2s+e8G2AEARAAAAIAQAAACAAAbWy0aDva244bYAQJAAABACAAAAQAgABkL8h5e0FGkAAAEAAAAgBAAAAIQPaCpr0gN8AIAuAXAAAgAAAEAIAAAFD2grIXJCMIgAAAEAAAAgBAAAAoe0HZC5IRBEAAAAgAAAEAIAAA1Pm3ILkBAAQAgABc1Q9n1E2/4fzgrQAAAABJRU5ErkJggg==)

```mathematica
In[14]:= DerivativeFilter[disk, {0, 2}]
Out[14]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABVElEQVR42u3cwQ2EIBBA0WGzhWllamVqZW4DciDBRfT9IzEB+ZnJBEdSRByBZnxsAQEEgAACQAABIIAAEEAA/se31cTDMJyOb9v2yHlFgBQEAggAAU+tgkqri16qoKvXKQKkIAJAAAHovQqapqmoWmh19lI6b+l7iQApCAQQAAKeWgXlzkb2fb+0GmlVBeXeq9YZkQiQgggAAQSgESkK/xFb1/V0fBzHV21crX0QAVIQASCAAPRSBR3H+eMppVdtXK19EAFSEAEggAAQQAAIIAAEEIC4cV9Qq36eu6E7WgoCAQSAgHjZF7FcV3BufJ7nrjcot/5a/76JACmIABBAAHqpgqJSt3Av9wVd3Q0uAqQgAkAAAYjO7wtaluWRd0eXvpcIkIJAAAEg4G1VUGn10st9QVevUwRIQQSAAAIQnX8RgwggAAQQAAIIAAEEgAACQMDt+QF6GGkJnPhDugAAAABJRU5ErkJggg==)

```mathematica
In[15]:= DerivativeFilter[zone, {0, 1}]
Out[15]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAFhElEQVR42u2d11IbQRBFRyByzhkKEG/8/69AUUWRC5FzDvKTx8dV02hWuyutrdtPXUI2u2v32ds9PT0l51zN1bGZmRnvLy4uer+np8f7r6+v3r+9vfX+zc2N9+/v773/9fXlsrByuRz0Ozo6gt//+PgI+lnZ5OSk9zc3N70/Pj7u/a2trT/X6WQtNf0DtNjK1g8YMvPz80HsvLy8eP/6+tr7l5eX3n96esr1Bkqlkvc7OzuDOOLnvb293v/8/AzeSxo08d53d3e9393dHbxmRYAQJAR5m5iY8P7Gxob3+/r6vH93d+f9arXq/dPT02L9z4IK6urqCvr8zvDwsPefn5+DaE1qR0dHQQRRASoChKD2ttL09LRPxCqViv/B4OBgMLHa39/3/vn5ectvwFI+9KncqILoE0004uji4iKYeCY1Il0RIAS1uQoaHR0NJl/EDt/mRcAOjYqCPpMs+lYNiigjIugTWVR9j4+Pia6ZSZ8iQAhqcwRRLby/v3v/5OSkaUmWlRzVarW6qEmKJt4jaz7f39/Ba6AaHBsbC6qsg4ODoGpSBAhBsroIYhgy0djb28v8l7EewuTIwg7Vi4WjpMa/h7iwfhdRw3oRFSO/v729rQgQgmTxCHp4eHChBfSsjCrirxpIBHaoTLiKlIdRHbEEzQSNCB0ZGXGhhfjl5WXvHx4eKgKEINnPCOIKF2sUaWxoaCjoU1EQO3n056QxXg9xxFpQf39/0GcPFRforQRNESAEtTmCssIOw5AKgQkXlYyFHaLJ8ptpXPmiSmRSZt373NycC/UIKQKEIFnd1sSkZoUkExmr3sKEy/JbhSAaVwlZmuaq4sDAgAvVi4QgIUiWOYJY5yF2WDOxlIy1YmX5RUAQr4E1NColPhMmoWz7vLq6UgQIQbL0CCJ2WCexlA9VjdXDw8+LpoJo3HjCZJbXyWdClciEThEgBAlBDRvrPDFbRGMSLgs7RUMQsfP29hZEKNWgpRIVAUKQEJTrznQLI/9K/ScmKbOUHp8JsUOVqAgQgoQg1+jOdIYYP7d6eNKomiLjKAab1rNSBAhBQlDDKEiqWGIwFfNni4YjCy8xmFIECEFCUCKztn8yxJhoWNhh2Fp+kZMy3ouVkPL62X2tYR1CkCyTcjTLsMSRpQRisBOT3BUBR9xswrK8tRrI7mjiSBEgBAlBDRvDir0xVmman1tzforcF0RjCyJxRFRy1YxDa/8qWev/oBAkBDVqnJPDljxraKo1ZDXGZ9gWIfli2yH7f6z2RbYjKgKEIFnmGzT4lrcmFqapC+W9Uz7GuBGDrYacd8TkNGboqyJACGpzBDGJSLNllUkZ1UvMTvmYVbNWIYgKh/vCmIjROHkgZtSbIkAIanMEcWc3w4qDJpIaExCiJmZeUBGMqoZbTfms+B2i++zsLNEERUWAENTmCKJKYT0nDYKcUS/iSlDM1MRmLsqze5lqhz4VEYeNsM4TM6ZMESAEyTyC+Dafmpry/urqqst6iCsRRD9mdnRW5WgmdFRlrO1Q7RA7XKHj4I6k2FEECEEyjyD2rlAJLCwsuNCCex6j7POYmmgt+hMpTDzpU6EROyy58zkkHVmvCBCCZD46+TbnCg6TsqWlpWBI/k/niFF9cVVL54gJQbJcEUSMEC88TZUKgWfKM/xbdZoq0UHsUNFZiR4VIBPDrE5Tpa2vr2uDhhAk+7kviGVV1mFWVlaCb/DZ2dlgyDfzTHkXsXudiR79PM6Up1E9EkHHx8eKACFI9nNrIt/+TDoYVlRHbNuj6mCfDFfHshqbH4OaZp7WwUN81tbWgupLp6kKQbLfVnLO1Rod1sqSNXHEkGe4URHRJ+KsE06LPC+Iz6FSqbjQKlu1WvX+zs6OIkAIkjnnnPsFScAeERegSawAAAAASUVORK5CYII=)

```mathematica
In[16]:= DerivativeFilter[rgb, {0, 1}]
Out[16]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAp0lEQVR42u3RIRIAIAwDwYDq/1+Lgw9UdgazJyOzK5WbroomO/28PfM3AAAACAAAAQAgAAAEAIAAABAAAAIAQAAACAAAAQAgAAAEAIAAABAAAAIAQAAACAAAAQAgAAAEAIAAABAAAAIAQAAACAAAAQAgAAAEAIAAABAAAAIAAIALAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAoJkeF7oCx5OTca4AAAAASUVORK5CYII=)

```mathematica
In[17]:= DerivativeFilter[vol, {0, 0, 1}]
Out[17]= -Image-
```

![12x10x8 volume](data:image/png;base64,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)

```mathematica
In[18]:= ImageChannels[DerivativeFilter[rgb, {0, 1}]]
Out[18]= 3

In[19]:= ImageDimensions[DerivativeFilter[vol, {0, 1, 0}]]
Out[19]= {12, 10, 8}
```

### Applications (2)

```mathematica
In[20]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];
```

```mathematica
In[21]:= Binarize[DerivativeFilter[zone, {0, 1}]]
Out[21]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABcUlEQVR42u3dQQ7DIAxEUVP1/ld2t9lESlSCCbxZV2rJaL5soHGLiIwbyrz18Witxa668qw+QaViQLG+sCMBDKA6tbMqaFXsXFnXE2s5+14JgCAIyup4zoagkWuXAAjaHEF5M5MrNVkz4EgCIAiC0t5OXRMqARAEQQk7dTiSAAiyF5QeQ12zJgEQFA7lVzyB6sbow2/7B0cSAEE0rAra4R5RrzVKAAQxgBjAAHp7I/ZEkzJ1+dipQZMACIKg8gjHZtvax7VLAARBkMpHAhhADGAAMYABxAAGEAMYQAxgADGAARQTb0c/fXNYAogBEBRjD6ljg3tEEgBB1B1Bu1VEvdYoARDEAGIAA2ilRmzViuiJv99KAARBUPgjRt0NcAmAoM0RdBYxb00cg1wJgCBVkBMr746GIHrLKEMTNEzQgCDqiCCH6eaIQRCFaaqmqRIEQZCZ8sNxJwEQBEHp1fR1jZ4EQBAE5Sr7Km88QZMACNpbPyXulrZcdWb6AAAAAElFTkSuQmCC)

### Properties & Relations (4)

```mathematica
In[22]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[23]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[24]:= ImageDimensions[DerivativeFilter[chk, {0, 1}]] === ImageDimensions[chk]
Out[24]= True

In[25]:= ImageChannels[DerivativeFilter[rgb, {1, 0}]] === ImageChannels[rgb]
Out[25]= True
```

### Neat Examples (2)

```mathematica
In[26]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];
```

```mathematica
In[27]:= DerivativeFilter[zone, {1, 1}]
Out[27]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAEkElEQVR42u2dW0/jMBCFXWgpdyQkHhD//6chQKAioNwv3Zfd7LeSj2I3dtJVzjyNoqhJ054vM2N7PJnP56vw246Pj/+4YWdnp/E/Pz8bf7lcNv7Ly0sY2iaTSeNvbW21+tvb21Gf56xWzSMJX19fjf/6+ho9J9dms9nf6wbboOYfYGCbEjsXFxfRkxaLRePf39//d1+SqKH85/N5FLk8//v7u/Gfnp4a/+HhYe37IdKtACNo5Aii9Gh3d3eNf3l5WeTtnxLJKJ/X/fn5iR4nLuhT8tPpNHp8b28v6vP5nJ6eRo/f3t5aAUaQbT0EUYaMdmpgR0UjKaghUroYEyv6Hx8fUTQdHh42/u7uboglrUzibm5urAAjyJaOINZ2mGR1wQ6RQtkyAlHYIWpqRFwpaHp8fIzeD5G4v78fxRTPSYmOrAAjaOQIKlVSJl6YyPA4jVJViRWPD2XPz8+tSFQ4YmSlakdWgBE0cgSVwg5lyCSLpmo1yleYGsqIazUSx6jv6OgoxEbTiCYrwAgygtZOshjtpGBH1WEUghh19JmU5UZHanCfx4kjRkRWgBFkBGWZqu2kYIdlXh5XiVguEml9IovREVGsEO15QUaQLQtBfJuXwk6pEa6Uwf0galChQimbCVfKNEgrwAgygkLObF5KWw2apyRcocdZ0//84yD/Gjhinef9/T3E5hFZAUaQrRVBuVMEU0rKm2CqjFzqPlUEqK5lBRhBRlCrVFPm7aSUlPs0XldFRLXNCDKCbEUQlDJvp0tJuTZ2Uu4n9/wU4zNR5XcrwAgat01CCKu2mklKFLEJCCoV7ZS6f94Py9FGkBFka0WQzQrwD2DrIRGrPa9m0yKTUuXrlPtMKXFbAUbQyBGkZFV73k7t2k6fbdBUWd4KMIJsgyNIYSQl+uoTOyk1HPqqXM+ysxVgBNmyEESJcUB5qMRn07CjENRlLZsVYASNHEGqjRgXYuS+2UPCqFnuiFv1f2JCf2n6CjtGkBFkWxtB7KWs1naVMoWjMFBPaSJX+USlmgVtBRhBtrUTMa5j4ipvrnuqjaMao2mqvEy8qAhQzecptfbNCjCCRo4gJg6UGxHEtzybmoYeZzt3WSmvkikV7fCcINZ/vb29WQFGkK07gtTbnBERm5HyHLbtChu27DSlnpOCHeK31D5ivK4VYASNHEGUFffJ4oY1bFOmFmIMta1h7paFuQPrfD5dklB2TXS/ICPIFk3E2FCUURA3rGGPaBWNDBUdqVE2NYCu+hoxySq1a+rZ2VkU11aAEeSV8lHjBjSUNpMy4khFIJRbjVJ2SFgcoQbNmWSVmopJ3J2fnzf+yclJiLUyswKMICOo1bg9H6WttvkjjhgJqF0kuqBJLY7gZ9Zug0/scFNsJrN8bkaQEWQrslKeb3bWOtSe7KqpKaMR+pvWj0glWYx2iB3a1dVVPML0f9AIMoJKfBBrR8SR6jut6jBqSqTy+zR+L9Z2iGKi8vr6unWjZyvACBo5gmaz2ar0ujA15U/NglbRjupN3Wc/an6Xg4ODEJs3xaiP0d1isWi9ZyvACBq3/QIfSe2Es1M6qwAAAABJRU5ErkJggg==)

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [ImageConvolve](../../image-processing/ImageConvolve/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
