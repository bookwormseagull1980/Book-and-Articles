"""Validate generated pages and their internal links before deployment."""
import json
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "_site"
PREFIX = "/Book-and-Articles/"
HOST = "bookwormseagull1980.github.io"


class PageParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self.ids = set()
        self.main = 0
        self.h1 = 0
        self.schema = False

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "main":
            self.main += 1
        if tag == "h1":
            self.h1 += 1
        if "id" in attrs:
            if attrs["id"] in self.ids:
                raise ValueError("Duplicate element ID: " + attrs["id"])
            self.ids.add(attrs["id"])
        if tag == "script" and attrs.get("type") == "application/ld+json":
            self.schema = True
        for attr in ("href", "src"):
            if attr in attrs:
                self.links.append(attrs[attr])


def check():
    pages = {}
    for path in SITE.rglob("*.html"):
        parser = PageParser()
        text = path.read_text(encoding="utf-8")
        parser.feed(text)
        assert parser.main == 1 and parser.h1 == 1, path
        assert parser.schema and 'rel="canonical"' in text, path
        assert "C:\\Users" not in text and "D:\\论文" not in text, path
        pages[path.relative_to(SITE).as_posix()] = parser
    errors = []
    for route, parser in pages.items():
        for link in parser.links:
            parts = urlsplit(link)
            if parts.netloc and parts.netloc != HOST:
                continue
            if parts.scheme in ("mailto", "javascript", "data"):
                continue
            if parts.path.startswith(PREFIX):
                target = unquote(parts.path[len(PREFIX):]) or "index.html"
            elif not parts.path:
                target = route
            else:
                continue
            if not (SITE / target).is_file():
                errors.append((route, link, "missing target"))
            elif parts.fragment and target in pages and unquote(parts.fragment) not in pages[target].ids:
                errors.append((route, link, "missing anchor"))
    catalogue = json.loads((SITE / "catalogue.json").read_text(encoding="utf-8"))
    assert len(catalogue) == 9
    assert len({a["doi"] for a in catalogue}) == 9
    assert len(pages) == 19, len(pages)
    home = (SITE / "index.html").read_text(encoding="utf-8")
    book = (SITE / "books/beyond-apparent-paths.html").read_text(encoding="utf-8")
    assert home.index('id="book"') < home.index('id="research"')
    assert '"@type": "Book"' in book
    assert "The Finger Pointing to the Moon" in book
    assert "Knowing Where to Stop Brings Stability" in book
    assert "five-page public preview" in book
    assert "Beyond%20Apparent%20Paths%20-%20Preview.pdf" in book
    for path in SITE.rglob("*.html"):
        text = path.read_text(encoding="utf-8")
        assert "Beyond%20Apparent%20Paths%20-%20The%20Practice" not in text, path
    if errors:
        raise AssertionError(errors)
    print(f"Validated {len(pages)} pages, internal links and 9 distinct DOI records.")


if __name__ == "__main__":
    check()
