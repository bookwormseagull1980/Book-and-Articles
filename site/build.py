"""Build the research portal from the collection's existing public metadata."""
from __future__ import annotations

import argparse
import html
import json
import posixpath
import re
import shutil
import subprocess
from datetime import date
from pathlib import Path
from urllib.parse import quote, unquote, urlsplit

import yaml
from markdown_it import MarkdownIt

ROOT = Path(__file__).resolve().parents[1]
GITHUB = "https://github.com/bookwormseagull1980/Book-and-Articles"
LEAN = "https://github.com/bookwormseagull1980/Lean-Formalization-of-Arithmetic-and-Basic-Physical-Interfaces"
V5 = "https://github.com/bookwormseagull1980/coarse-graining-genesis"
DATA = json.loads((ROOT / "site/research.json").read_text(encoding="utf-8"))
ROUTES = {"README.md": "index.html", "CITATION.md": "citation.html"}
for topic in DATA["topics"]:
    ROUTES[f"topics/{topic['slug']}.md"] = f"topics/{topic['slug']}.html"
for slug in DATA["reading_order"]:
    ROUTES[f"articles/{slug}/README.md"] = f"papers/{slug}.html"
ROUTES["books/beyond-apparent-paths/README.md"] = "books/beyond-apparent-paths.html"


def escape(value):
    return html.escape(str(value), quote=True)


def rewrite_link(value, source, base):
    parts = urlsplit(value)
    if parts.scheme or parts.netloc or not parts.path:
        return value
    target = posixpath.normpath(posixpath.join(posixpath.dirname(source), unquote(parts.path)))
    suffix = ("?" + parts.query if parts.query else "") + ("#" + parts.fragment if parts.fragment else "")
    if target in ROUTES:
        return base + ROUTES[target] + suffix
    if target.startswith("site/assets/"):
        return base + "assets/" + quote(target.removeprefix("site/assets/")) + suffix
    return GITHUB + "/blob/main/" + quote(target, safe="/") + suffix


def render_markdown(text, source, base, drop_title=False):
    parser = MarkdownIt("commonmark", {"html": False, "linkify": False}).enable("table")
    env = {}
    tokens = parser.parse(text, env)
    seen = {}
    for token in tokens:
        if token.type == "heading_open":
            index = tokens.index(token)
            inline = tokens[index + 1]
            anchor = re.sub(r"[^\w\s-]", "", inline.content.lower())
            anchor = re.sub(r"[\s_]+", "-", anchor).strip("-")
            count = seen.get(anchor, 0)
            seen[anchor] = count + 1
            token.attrSet("id", anchor + (f"-{count}" if count else ""))
        for child in token.children or []:
            if child.type == "link_open":
                child.attrSet("href", rewrite_link(child.attrGet("href"), source, base))
            elif child.type == "image":
                child.attrSet("src", rewrite_link(child.attrGet("src"), source, base))
    if drop_title and tokens and tokens[0].type == "heading_open" and tokens[0].tag == "h1":
        tokens = tokens[3:]
    return parser.renderer.render(tokens, parser.options, env)


def text_excerpt(text, limit=220):
    text = re.sub(r"\s+", " ", text).strip()
    return text if len(text) <= limit else text[:limit].rsplit(" ", 1)[0] + "..."


def publication_links(identifiers):
    links = []
    for label, doi_key, url_key in [
        ("arXiv", "arxiv_doi", "arxiv_url"),
        ("SSRN", "ssrn_doi", "ssrn_url"),
    ]:
        doi = identifiers.get(doi_key)
        url = identifiers.get(url_key)
        if doi or url:
            if not url:
                url = ("https://arxiv.org/abs/" + doi.split("arXiv.", 1)[1]
                       if label == "arXiv" else "https://ssrn.com/abstract=" + doi.rsplit(".", 1)[1])
            links.append({"label": label, "url": url, "doi": doi})
    doi = identifiers.get("zenodo_concept_doi") or identifiers["zenodo_doi"]
    links.append({"label": "Zenodo", "url": "https://doi.org/" + doi, "doi": doi})
    return links


def paper_links(article):
    return "\n".join(
        f'<a href="{escape(link["url"])}">{escape(link["label"])}</a>'
        for link in article["publications"]
    )


def main():
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("--output", default="_site")
    args_parser.add_argument("--site-url", default="https://bookwormseagull1980.github.io/Book-and-Articles/")
    args = args_parser.parse_args()
    output = (ROOT / args.output).resolve()
    if ROOT not in output.parents:
        raise ValueError("The output directory must be inside the repository.")
    output.mkdir(parents=True, exist_ok=True)
    base = args.site_url.rstrip("/") + "/"
    paths = []
    articles = []
    tree = set(subprocess.check_output(["git", "-c", f"safe.directory={ROOT}", "ls-tree", "-r", "--name-only", "HEAD"], cwd=ROOT, text=True).splitlines())
    for slug in DATA["reading_order"]:
        directory = ROOT / "articles" / slug
        identifiers = yaml.safe_load((directory / "metadata/IDENTIFIERS.yml").read_text(encoding="utf-8"))
        citation = yaml.safe_load((directory / "metadata/CITATION.cff").read_text(encoding="utf-8"))
        pdf_path = f"articles/{slug}/{identifiers['pdf']['filename']}"
        if pdf_path not in tree:
            raise ValueError(f"Untracked PDF reference: {pdf_path}")
        publications = publication_links(identifiers["identifiers"])
        articles.append({
            "slug": slug, "title": identifiers["title"],
            "doi": publications[0]["doi"],
            "archive_doi": publications[-1]["doi"],
            "primary_url": publications[0]["url"],
            "publications": publications,
            "ssrn": identifiers["identifiers"].get("ssrn_doi"),
            "arxiv": identifiers["identifiers"].get("arxiv_doi"),
            "pdf": "https://raw.githubusercontent.com/bookwormseagull1980/Book-and-Articles/main/" + quote(pdf_path, safe="/"),
            "github": GITHUB + "/tree/main/articles/" + slug,
            "keywords": citation.get("keywords", []),
            "abstract": citation.get("abstract", ""),
            "date": str(citation.get("date-released", "")),
            "markdown": (directory / "README.md").read_text(encoding="utf-8"),
        })

    def page(route, title, content, description=None, article=None):
        description = description or DATA["description"]
        canonical = base + route
        meta = ""
        structured = {
            "@context": "https://schema.org", "@type": "WebSite",
            "name": DATA["title"], "url": base,
            "author": {"@type": "Person", "name": "Jinku Guo", "identifier": "https://orcid.org/0009-0000-6600-6171"},
        }
        if article:
            for name, value in [("citation_title", article["title"]), ("citation_author", "Jinku Guo"),
                                ("citation_doi", article["doi"]), ("citation_pdf_url", article["pdf"]),
                                ("citation_publication_date", article["date"])]:
                if value:
                    meta += f'<meta name="{name}" content="{escape(value)}">\n'
            structured = {
                "@context": "https://schema.org", "@type": "ScholarlyArticle",
                "headline": article["title"], "url": canonical,
                "identifier": "https://doi.org/" + article["doi"],
                "sameAs": [link["url"] for link in article["publications"]],
                "author": {"@type": "Person", "name": "Jinku Guo"},
                "keywords": article["keywords"], "abstract": article["abstract"],
                "encoding": {"@type": "MediaObject", "contentUrl": article["pdf"], "encodingFormat": "application/pdf"},
            }
            if re.fullmatch(r"\d{4}-\d{2}-\d{2}", article["date"]):
                structured["datePublished"] = article["date"]
        document = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{escape(title)} | Jinku Guo</title>
<meta name="description" content="{escape(description)}">
<link rel="canonical" href="{escape(canonical)}">
<meta property="og:type" content="{'article' if article else 'website'}">
<meta property="og:title" content="{escape(title)}">
<meta property="og:description" content="{escape(description)}">
<meta property="og:url" content="{escape(canonical)}">
<meta property="og:image" content="{base}assets/book-and-articles-social.png">
<meta name="twitter:card" content="summary_large_image">
{meta}<script type="application/ld+json">{json.dumps(structured, ensure_ascii=False).replace("</", "<\\/")}</script>
<link rel="stylesheet" href="{base}assets/style.css">
<script defer src="{base}assets/search.js"></script>
<script>window.MathJax={{tex:{{inlineMath:[['$','$'],['\\\\(','\\\\)']]}}}};</script>
<script async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
</head><body>
<a class="skip" href="#main">Skip to content</a>
<header class="site-header"><div class="container header-inner">
<a class="brand" href="{base}index.html">Foundations from Minimal Change</a>
<nav aria-label="Main navigation"><a href="{base}index.html#research">Research</a>
<a href="{base}index.html#papers">Papers</a><a href="{base}verification.html">Proofs &amp; code</a>
<a href="{base}citation.html">Cite</a><a href="{base}author.html">Author</a></nav></div></header>
<main id="main" class="container">{content}</main>
<footer class="container"><p>Jinku Guo · Northwestern Polytechnical University ·
<a href="https://orcid.org/0009-0000-6600-6171">ORCID</a> ·
<a href="{GITHUB}">GitHub collection</a> ·
<a href="{GITHUB}/discussions">Research discussions</a></p>
<p>Preprints are available on arXiv and SSRN. Zenodo supplies the archival records.</p></footer>
</body></html>'''
        destination = output / route
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(document, encoding="utf-8")
        paths.append(route)

    topic_rows = "".join(
        f'<li><a href="{base}topics/{t["slug"]}.html">{escape(t["short"])}</a></li>'
        for t in DATA["topics"]
    )
    article_rows = []
    for number, article in enumerate(articles, 1):
        search = escape(" ".join([article["title"], article["abstract"], *article["keywords"]]))
        keywords = escape(" · ".join(article["keywords"][:6]))
        article_rows.append(f'''<article class="paper-row" data-search="{search}">
<span class="paper-index">{number:02d}</span><div>
<h3><a href="{base}papers/{article["slug"]}.html">{escape(article["title"])}</a></h3>
<p>{escape(text_excerpt(article["abstract"]))}</p>
<p class="paper-keywords">{keywords}</p>
<p class="paper-links">{paper_links(article)}
<a href="{article["pdf"]}">Read PDF</a>
<a href="{article["github"]}">Files &amp; supplement</a></p></div></article>''')

    landing = f'''<section class="opening"><p class="eyebrow">Research by Jinku Guo</p>
<h1>Foundations from<br>Minimal Change</h1><p class="subtitle">Mathematics, physics and computation</p>
<p class="lead">A common generating source connects arithmetic histories,
compact gauge structure, long-range response and finite certified computation.</p>
<p>The programme follows what is produced by a construction and what must be
retained when it is observed at finite resolution. Its results include a
history-dependent natural-number system, quantitative prime projectors, a generated
Yang-Mills spectral-gap theorem, and polynomial-cost trajectory certificates.</p>
<p class="opening-links"><a href="#papers">Browse the papers <span aria-hidden="true">↓</span></a>
<a href="{base}verification.html">Explore proofs and calculations <span aria-hidden="true">→</span></a></p>
</section>
<section id="research" class="research-section"><h2>A common source, several research directions</h2>
<figure><img id="research-map" src="{base}assets/research-map.png" usemap="#research-navigation" alt="Minimal change leads to Gaussian correlation and retained spectral data, then branches into arithmetic histories, physical response and certified computation." width="1600" height="720">
<map name="research-navigation">
<area shape="rect" coords="475,20,1125,112" data-coords="475,20,1125,112" href="{base}topics/generated-arithmetic.html" alt="Minimal change and the arithmetic source" title="Minimal change and the arithmetic source">
<area shape="rect" coords="475,158,1125,254" data-coords="475,158,1125,254" href="{base}topics/generated-arithmetic.html" alt="Gaussian correlation and the generated carrier" title="Gaussian correlation and the generated carrier">
<area shape="rect" coords="475,300,1125,396" data-coords="475,300,1125,396" href="{base}topics/prime-operators.html" alt="Retained spectral data and prime operators" title="Retained spectral data and prime operators">
<area shape="rect" coords="20,495,510,552" data-coords="20,495,510,552" href="{base}topics/generated-arithmetic.html" alt="Arithmetic and histories" title="Arithmetic and histories">
<area shape="rect" coords="20,552,510,615" data-coords="20,552,510,615" href="{base}topics/prime-operators.html" alt="Prime directions and Hodge operators" title="Prime directions and Hodge operators">
<area shape="rect" coords="20,615,510,680" data-coords="20,615,510,680" href="{base}topics/generated-arithmetic.html" alt="Value response and memory" title="Value response and memory">
<area shape="rect" coords="555,495,1045,615" data-coords="555,495,1045,615" href="{base}topics/gauge-and-mass-gap.html" alt="Gauge structure and the generated Yang-Mills spectral gap" title="Gauge structure and the generated Yang-Mills spectral gap">
<area shape="rect" coords="555,615,1045,680" data-coords="555,615,1045,680" href="{base}topics/gravity-and-cosmology.html" alt="Gravity and cosmic structure" title="Gravity and cosmic structure">
<area shape="rect" coords="1090,495,1580,680" data-coords="1090,495,1580,680" href="{base}topics/certified-computation.html" alt="Certified computation and polynomial bit cost" title="Certified computation and polynomial bit cost">
</map>
<figcaption>The branches share generating data. Their theorems specify the maps,
realizations and response classes connecting one stage to the next.
<a href="{base}assets/research-map.png">Full-resolution diagram</a>.</figcaption></figure>
<ul class="topic-list">{topic_rows}</ul></section>
<section class="reading"><h2>Reading routes</h2><div class="reading-columns">
<div><h3>Mathematical construction</h3><p>Begin with natural-number arithmetic,
then prime decomposition and generated spectra. Expression histories connect
arithmetic values to later response and memory.</p></div>
<div><h3>Physical response</h3><p>Read Compact Internal Space I and II, then
the gravity criterion and coarse-graining. The third compact-space paper develops
finite-window cosmic structure.</p></div>
<div><h3>Formal and numerical evidence</h3><p>The Lean library gives declarations
and theorem-to-paper mappings. V5 gives executable source records, interval
certificates and observational comparisons.</p></div></div></section>
<section id="papers"><div class="section-heading"><h2>Papers</h2>
<label class="search-label">Find a paper<input id="paper-search" type="search" placeholder="Title, subject or result" aria-controls="paper-list"></label></div>
<p id="search-status" class="muted" aria-live="polite">9 papers</p>
<div id="paper-list">{"".join(article_rows)}</div></section>
<section><h2>Books</h2><p><a href="{base}books/beyond-apparent-paths.html">Beyond Apparent Paths:
The Practice of Seeking Truth from Facts</a></p></section>'''
    page("index.html", DATA["title"] + ": " + DATA["subtitle"], landing)

    for topic in DATA["topics"]:
        source = f"topics/{topic['slug']}.md"
        markdown = (ROOT / source).read_text(encoding="utf-8")
        title = markdown.splitlines()[0].removeprefix("# ")
        page(f"topics/{topic['slug']}.html", title,
             '<article class="prose">' + render_markdown(markdown, source, base) + "</article>")

    for article in articles:
        tools = f'''<p class="article-actions">{paper_links(article)}
<a href="{article["pdf"]}">Read manuscript PDF</a>
<a href="{article["github"]}">Files and citation</a></p>'''
        page(f"papers/{article['slug']}.html", article["title"],
             '<article class="prose paper-detail"><p class="eyebrow">Research article · Jinku Guo</p><h1>'
             + escape(article["title"]) + "</h1>" + tools
             + render_markdown(article["markdown"], f"articles/{article['slug']}/README.md", base, drop_title=True)
             + "</article>", text_excerpt(article["abstract"], 250), article)
    for source, route in [("CITATION.md", "citation.html"), ("books/beyond-apparent-paths/README.md", "books/beyond-apparent-paths.html"),
                          ("site/verification.md", "verification.html"), ("site/author.md", "author.html")]:
        markdown = (ROOT / source).read_text(encoding="utf-8")
        page(route, markdown.splitlines()[0].removeprefix("# "),
             '<article class="prose">' + render_markdown(markdown, source, base) + "</article>")

    shutil.copytree(ROOT / "site/assets", output / "assets", dirs_exist_ok=True)
    (output / ".nojekyll").write_text("", encoding="utf-8")
    (output / "robots.txt").write_text(f"User-agent: *\nAllow: /\nSitemap: {base}sitemap.xml\n", encoding="utf-8")
    entries = "".join(f"<url><loc>{escape(base + path)}</loc></url>" for path in paths)
    (output / "sitemap.xml").write_text('<?xml version="1.0" encoding="UTF-8"?><urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">' + entries + "</urlset>", encoding="utf-8")
    (output / "catalogue.json").write_text(json.dumps([{k: v for k, v in a.items() if k != "markdown"} for a in articles], ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Built {len(paths)} pages, {len(articles)} paper records and sitemap in {output}")


if __name__ == "__main__":
    main()
