# stb (vendored, unmodified)

`stb_image.h` v2.30 and `stb_image_write.h` v1.16 by Sean Barrett — public domain
(alternatively MIT), <http://nothings.org/stb>.

Vendored rather than taken as a system dependency because a decoder is not optional
for `Import["photo.jpg"]`: a build without one turns the whole image subsystem into
something that can only read arrays typed by hand. Two headers with no dependencies
beyond libc are cheaper than making libpng and libjpeg hard requirements, and cheaper
still than the alternative of writing a baseline JPEG decoder.

**Do not edit these files.** Warnings are suppressed at the single point of inclusion
(`src/imageio.c`) so that the upstream text stays byte-identical and can be re-fetched:

    curl -O https://raw.githubusercontent.com/nothings/stb/master/stb_image.h
    curl -O https://raw.githubusercontent.com/nothings/stb/master/stb_image_write.h
