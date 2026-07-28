#!/usr/bin/env python3
"""Rasterise Mathilda's DensityPlot scenes (Plotly JSON from the REPL's pipe
mode) into one composite PNG.  The polygons, their coordinates and their
colours are all Mathilda's; this only draws them."""
import sys, json, re, zlib, struct

PANEL   = 200
GAP     = 16
MARGIN_L, MARGIN_R = 26, 26
MARGIN_T, MARGIN_B = 24, 52
BG      = (255, 255, 255)
INK     = (60, 60, 60)

FONT = {  # 5x7, one string of 5 chars per row, '#' = ink
 '0':["01110","10001","10011","10101","11001","10001","01110"],
 '1':["00100","01100","00100","00100","00100","00100","01110"],
 '2':["01110","10001","00001","00010","00100","01000","11111"],
 '3':["11111","00010","00100","00010","00001","10001","01110"],
 '4':["00010","00110","01010","10010","11111","00010","00010"],
 '5':["11111","10000","11110","00001","00001","10001","01110"],
 '6':["00110","01000","10000","11110","10001","10001","01110"],
 '7':["11111","00001","00010","00100","01000","01000","01000"],
 '8':["01110","10001","10001","01110","10001","10001","01110"],
 '9':["01110","10001","10001","01111","00001","00010","01100"],
 '.':["00000","00000","00000","00000","00000","01100","01100"],
 '-':["00000","00000","00000","11111","00000","00000","00000"],
 '+':["00000","00100","00100","11111","00100","00100","00000"],
 '=':["00000","00000","11111","00000","11111","00000","00000"],
 't':["01000","01000","11110","01000","01000","01001","00110"],
 'u':["00000","00000","10001","10001","10001","10011","01101"],
 'x':["00000","00000","10001","01010","00100","01010","10001"],
 'y':["00000","00000","10001","10001","10011","01101","00001"],
 ' ':["00000"]*7,
}

class Img:
    def __init__(self, w, h, bg):
        self.w, self.h = w, h
        self.px = bytearray(bg * (w * h))
    def rect(self, x0, y0, x1, y1, c):
        x0 = max(0, int(x0)); y0 = max(0, int(y0))
        x1 = min(self.w, int(x1)); y1 = min(self.h, int(y1))
        row = bytes(c) * max(0, x1 - x0)
        for y in range(y0, y1):
            o = (y * self.w + x0) * 3
            self.px[o:o + len(row)] = row
    def frame(self, x0, y0, x1, y1, c, t=1):
        self.rect(x0, y0, x1, y0 + t, c); self.rect(x0, y1 - t, x1, y1, c)
        self.rect(x0, y0, x0 + t, y1, c); self.rect(x1 - t, y0, x1, y1, c)
    def text(self, x, y, s, c, k=2):
        for ch in s:
            g = FONT.get(ch, FONT[' '])
            for r, rowbits in enumerate(g):
                for col, b in enumerate(rowbits):
                    if b == '1':
                        self.rect(x + col * k, y + r * k, x + (col + 1) * k, y + (r + 1) * k, c)
            x += 6 * k
        return x
    def width_of(self, s, k=2):
        return 6 * k * len(s)
    def png(self):
        raw = bytearray()
        for y in range(self.h):
            raw.append(0)
            raw += self.px[y * self.w * 3:(y + 1) * self.w * 3]
        def chunk(tag, data):
            return (struct.pack(">I", len(data)) + tag + data
                    + struct.pack(">I", zlib.crc32(tag + data) & 0xffffffff))
        return (b"\x89PNG\r\n\x1a\n"
                + chunk(b"IHDR", struct.pack(">IIBBBBB", self.w, self.h, 8, 2, 0, 0, 0))
                + chunk(b"IDAT", zlib.compress(bytes(raw), 9))
                + chunk(b"IEND", b""))

def scenes(path):
    """Every {"type":"plot"} record, in order."""
    out = []
    for line in open(path):
        line = line.strip()
        if '"type":"plot"' not in line:
            continue
        out.append(json.loads(line)["payload"]["data"])
    return out

RGBA = re.compile(r"rgba\((\d+),(\d+),(\d+)")

def draw_panel(img, ox, oy, traces):
    """Each trace is one axis-aligned cell of Mathilda's density grid."""
    for tr in traces:
        xs, ys = tr.get("x"), tr.get("y")
        if not xs or not ys:
            continue
        m = RGBA.match(tr.get("fillcolor") or "")
        if not m:
            continue
        c = tuple(int(v) for v in m.groups())
        x0, x1 = min(xs), max(xs)
        y0, y1 = min(ys), max(ys)
        # data [0,1]^2 -> panel pixels, y flipped (image rows run downward)
        px0 = ox + x0 * PANEL; px1 = ox + x1 * PANEL
        py0 = oy + (1 - y1) * PANEL; py1 = oy + (1 - y0) * PANEL
        img.rect(px0, py0, max(px1, px0 + 1), max(py1, py0 + 1), c)

def colour(z):
    """The same ramp the ColorFunction uses, for the legend bar."""
    return (int(255 * min(1.0, 1.0 + z)),
            int(255 * (1.0 - abs(z))),
            int(255 * min(1.0, 1.0 - z)))

def main(src, dst, labels):
    data = scenes(src)
    n = len(data)
    assert n == len(labels), "%d scenes, %d labels" % (n, len(labels))
    W = MARGIN_L + n * PANEL + (n - 1) * GAP + MARGIN_R
    H = MARGIN_T + PANEL + MARGIN_B
    img = Img(W, H, BG)
    for i, traces in enumerate(data):
        ox = MARGIN_L + i * (PANEL + GAP)
        draw_panel(img, ox, MARGIN_T, traces)
        img.frame(ox - 1, MARGIN_T - 1, ox + PANEL + 1, MARGIN_T + PANEL + 1, INK)
        lab = labels[i]
        img.text(ox + (PANEL - img.width_of(lab)) // 2, 5, lab, INK)

    # legend: the ramp from -1 to +1, centred under the panels
    bw, bh = 300, 11
    bx = (W - bw) // 2
    by = MARGIN_T + PANEL + 22
    for k in range(bw):
        img.rect(bx + k, by, bx + k + 1, by + bh, colour(-1.0 + 2.0 * k / (bw - 1)))
    img.frame(bx - 1, by - 1, bx + bw + 1, by + bh + 1, INK)
    img.text(bx - 30, by + 1, "-1", INK)
    img.text(bx + bw // 2 - 6, by + bh + 6, "0", INK)
    img.text(bx + bw + 8, by + 1, "+1", INK)
    img.text(bx - 62, by + 1, "u", INK)
    open(dst, "wb").write(img.png())
    print("%s  %dx%d  %d bytes" % (dst, W, H, len(open(dst, 'rb').read())))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
