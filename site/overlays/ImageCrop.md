### Worked examples

```mathematica
In[1]:= ImageDimensions[ImageCrop[Image[{{0., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}}], {1, 1}]]  (* crop to 1x1 about the centre *)
```

```mathematica
In[1]:= ImageData[ImageCrop[Image[{{0., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}}], {1, 1}]]  (* the centre pixel survives *)
```

```mathematica
In[1]:= ImageDimensions[ImageCrop[Image[{{0., 0., 0., 0.}, {0., 1., 1., 0.}, {0., 0., 0., 0.}}]]]  (* trim the uniform black border *)
```

### Notes

`ImageCrop[image, {w, h}]` crops to `w × h` about the centre, any odd remainder going to the
right and bottom — the same floor-division convention the kernel centres use, which is what
makes `ImageCrop[ImagePad[image, m], ImageDimensions[image]]` exactly the original image. A
crop may not enlarge.

`ImageCrop[image]` with no size trims a uniform border, asking how much of the frame carries no
information. The border colour is read from a corner rather than assumed black, so a scanned
page's white margin is trimmed too. An entirely uniform image comes back unchanged — there is
no content to keep and a zero-sized image is not an image.
