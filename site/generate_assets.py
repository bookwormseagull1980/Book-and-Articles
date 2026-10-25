"""Create the static research diagram and repository sharing images with Pillow."""
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

OUT = Path(__file__).resolve().parent / "assets"


def font(size, bold=False):
    names = [
        Path("C:/Windows/Fonts") / ("segoeuib.ttf" if bold else "segoeui.ttf"),
        Path("/usr/share/fonts/truetype/dejavu") / ("DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"),
    ]
    return ImageFont.truetype(str(next(p for p in names if p.is_file())), size)


def centre(draw, bounds, lines, size=25, color="#244b43"):
    x0, y0, x1, y1 = bounds
    height = len(lines) * (size + 10)
    for n, line in enumerate(lines):
        width = draw.textlength(line, font=font(size))
        draw.text(((x0 + x1 - width) / 2, (y0 + y1 - height) / 2 + n * (size + 10)), line, font=font(size), fill=color)


def arrow(draw, start, end, color="#788e84"):
    draw.line((start, end), fill=color, width=3)
    x, y = end
    draw.polygon(((x - 8, y - 11), (x + 8, y - 11), (x, y)), fill=color)


def diagram():
    im = Image.new("RGB", (1600, 720), "#ffffff")
    d = ImageDraw.Draw(im)
    boxes = [
        ((475, 20, 1125, 112), ["Minimal unbiased compositional change"], "#eaf3ed"),
        ((475, 158, 1125, 254), ["Gaussian correlation and finite observation"], "#eaf3ed"),
        ((475, 300, 1125, 396), ["Retained spectral data and response"], "#edf2f7"),
        ((20, 495, 510, 680), ["Arithmetic and histories", "Prime directions / Hodge operators", "Value response and memory"], "#edf2f7"),
        ((555, 495, 1045, 680), ["Gauge and matter structures", "Generated Yang-Mills spectral gap", "Gravity and cosmic structure"], "#f7eef0"),
        ((1090, 495, 1580, 680), ["Certified computation", "Finite records / residual bounds", "Polynomial bit cost"], "#eaf3ed"),
    ]
    for bounds, lines, fill in boxes:
        d.rounded_rectangle(bounds, radius=6, fill=fill, outline="#bacbc4", width=2)
        centre(d, bounds, lines, size=24 if len(lines) > 1 else 26)
    arrow(d, (800, 112), (800, 158))
    arrow(d, (800, 254), (800, 300))
    d.line(((800, 396), (800, 447)), fill="#788e84", width=3)
    d.line(((265, 447), (1335, 447)), fill="#788e84", width=3)
    for x in (265, 800, 1335):
        arrow(d, (x, 447), (x, 495))
    im.save(OUT / "research-map.png", optimize=True)


def social(filename, eyebrow, title, lines, tags, accent):
    im = Image.new("RGB", (1280, 640), "#ffffff")
    d = ImageDraw.Draw(im)
    d.rectangle((0, 0, 15, 640), fill=accent)
    d.text((70, 49), eyebrow, font=font(21, True), fill=accent)
    y = 109
    for line in title:
        d.text((70, y), line, font=font(51, True), fill="#25362f")
        y += 67
    y += 24
    for line in lines:
        d.text((73, y), line, font=font(26), fill="#42594f")
        y += 41
    d.line((73, 516, 1205, 516), fill="#cdd9d3", width=2)
    d.text((73, 548), tags, font=font(20), fill=accent)
    d.text((73, 585), "Jinku Guo  |  Mathematics, physics and computation", font=font(18), fill="#617169")
    im.save(OUT / filename, optimize=True)
    assert (OUT / filename).stat().st_size < 1_000_000


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    diagram()
    social("book-and-articles-social.png", "FOUNDATIONS FROM MINIMAL CHANGE",
           ["A common source.", "Mathematics, physics, computation."],
           ["Generated arithmetic and retained histories", "Gauge structure, spectral gaps and cosmic response", "Certified finite computation and Lean proofs"],
           "Arithmetic  /  Yang-Mills  /  Gravity  /  P vs NP", "#11685e")
    social("lean-social.png", "LEAN 4 / MATHLIB",
           ["Formal foundations", "from minimal change"],
           ["Generated arithmetic and history dynamics", "Gauge-matter structure and Yang-Mills operators", "Kernel-checked theorems and explicit axiom audits"],
           "Lean 4  /  Gaussian source  /  Arithmetic  /  Spectral gap", "#235a87")
    social("v5-social.png", "COARSE-GRAINING GENESIS FRAMEWORK",
           ["One generating source.", "Correlated physical responses."],
           ["Spectral geometry, gauge structure and mass scales", "Emergent gravity and cosmological readouts", "Reproducible records and interval certificates"],
           "Yang-Mills  /  Particle physics  /  Gravity  /  Cosmology", "#93394d")
    print("Generated research diagram and three sharing images.")
