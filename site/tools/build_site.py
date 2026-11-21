#!/usr/bin/env python3
"""Build the QMCkl static site from org-exported HTML.

Usage: build_site.py REPO_DIR OUT_DIR

Per documentation page (org/*.html, order = org/table_of_contents):
  1. mask <pre>/<code> regions (code is never treated as math)
  2. convert \\( ... \\) and \\[ ... \\] TeX to native MathML (latex2mathml)
  3. restore code regions, highlight with Pygments (classes styled by CSS)
  4. rewrite <head>, strip org-info.js scripts / postamble / home div
  5. build rustdoc-style layout: fixed sidebar (cross-page nav + page TOC)
  6. fix relative links and audit the result

Code visibility policy (user, 2026-11-21): the API pages show only the
qmckl.h C declarations (prototypes, typedefs, enums, defines) and all
Fortran kernels. C function bodies, usage fragments, shell/plain blocks
and the qmckl_compute_* C prototypes are hidden. The worked-example page
(qmckl_examples.html) keeps all of its C/Fortran/shell blocks — complete
programs — and links each block to the downloadable program file it
belongs to (OUT_DIR/code/): mo_ortho.c and ao_grid.f90 come from
org-babel-tangle (run by export_docs.sh, noweb-aware), the other two are
the verbatim source blocks of their (untangled) sections.

Also highlights <pre data-lang="..."> blocks in the hand-written landing
page (tools/index.html, next to this script) and writes OUT_DIR/index.html.

Build-time only: the published site contains no JavaScript at all.
"""
import html
import json
import re
import sys
from pathlib import Path

import latex2mathml.converter as l2m
from pygments import highlight
from pygments.formatters import HtmlFormatter
from pygments.lexers import get_lexer_by_name

if len(sys.argv) != 3:
    sys.exit(__doc__)
REPO = Path(sys.argv[1]).resolve()
SITE = Path(sys.argv[2]).resolve()  # output: receives docs/, code/, index.html
ORG = REPO / "org"
DOCS_OUT = SITE / "docs"
TOOLSDIR = Path(__file__).resolve().parent

TOC = [l.strip() for l in (ORG / "table_of_contents").read_text().splitlines() if l.strip()]

# org "src src-X" language -> (pygments lexer, display label)
LEXERS = {
    "c": ("c", "C"), "C": ("c", "C"),
    "f90": ("fortran", "Fortran"), "fortran": ("fortran", "Fortran"), "f": ("fortran", "Fortran"),
    "python": ("python", "Python"), "py": ("python", "Python"),
    "bash": ("bash", "Shell"), "sh": ("bash", "Shell"), "shell": ("bash", "Shell"),
    "emacs-lisp": ("emacs-lisp", "Emacs Lisp"), "elisp": ("emacs-lisp", "Emacs Lisp"),
    "makefile": ("makefile", "Makefile"),
}

PRE_RE = re.compile(r"<pre\b[^>]*>.*?</pre>", re.S)
CODE_RE = re.compile(r"<code\b[^>]*>.*?</code>", re.S)
# defensive: org fontifies src blocks with inline-CSS spans even when the
# export was told not to (emacs/org versions differ); strip face wrappers
SPAN_RE = re.compile(r"</?span\b[^>]*>")
MATH_D_RE = re.compile(r"\\\[.*?\\\]", re.S)
MATH_I_RE = re.compile(r"\\\((.*?)\\\)", re.S)
# raw LaTeX eqnarray blocks written in the org body (not as $…$ math):
# org exports them as plain text. eqnarray ≡ array{rcl} for the simple
# "X & = & Y \\ & = & Z" alignment these docs use, and latex2mathml
# renders array (real mtable) but not eqnarray (dumps literal '&')
EQA_RE = re.compile(r"\\begin\{eqnarray\*?\}.*?\\end\{eqnarray\*?\}", re.S)


# ---------------------------------------------------------------- helpers

def extract_balanced(text: str, start: int) -> tuple[int, int, int]:
    """For the <div opening at `start`, return (content_start, content_end, region_end):
    content spans (content_start, content_end) — exclusive of the closing tag;
    the whole div region (incl. the closing tag) spans (start, region_end)."""
    depth, i = 0, start
    while True:
        m = re.search(r"<div\b|</div>", text[i:], re.I)
        if not m:
            raise ValueError("unbalanced div")
        i += m.end()
        depth += -1 if m.group(0).startswith("</") else 1
        if depth == 0:
            return text.index(">", start) + 1, i - len("</div>"), i


def first_paragraph(html_text: str, limit: int = 155) -> str:
    m = re.search(r"<p>(.*?)</p>", html_text, re.S)
    if not m:
        return "QMCkl kernel library documentation."
    t = re.sub(r"<[^>]+>", "", m.group(1))
    t = re.sub(r"\s+", " ", html.unescape(t)).strip()
    if len(t) > limit:
        t = t[: limit - 1].rsplit(" ", 1)[0] + "…"
    return t


# ---------------------------------------------------------------- math

def render_math(tex: str, display: bool, audit: dict, name: str) -> str:
    try:
        mml = l2m.convert(tex)
    except Exception as e:  # keep visible, flag it
        audit["math_fail"].append({"file": name, "tex": " ".join(tex.split())[:90],
                                   "error": f"{type(e).__name__}: {e}"[:90]})
        raw = f'<code class="math-raw">{html.escape(tex)}</code>'
        return f'<div class="math-display">{raw}</div>' if display else raw
    if display:
        mml = mml.replace('display="inline"', 'display="block"', 1)
        return f'<div class="math-display" data-tex="{html.escape(tex, quote=True)}">{mml}</div>'
    return f'<span class="math-inline" data-tex="{html.escape(tex, quote=True)}">{mml}</span>'


def pass_math(text: str, audit: dict, name: str) -> str:
    # raw eqnarray blocks first: rewrite to array{rcl}, convert as display
    # math; they sit between <p>s as bare text, so the replacement <div>
    # is a valid block sibling
    def eqa_sub(m):
        tex = re.sub(r"\\begin\{eqnarray\*?\}", r"\\begin{array}{rcl}", m.group(0))
        tex = re.sub(r"\\end\{eqnarray\*?\}", r"\\end{array}", tex)
        audit["math_eqnarray"] += 1
        return render_math(tex, True, audit, name)

    text = EQA_RE.sub(eqa_sub, text)
    n_d, n_i = len(MATH_D_RE.findall(text)), len(MATH_I_RE.findall(text))
    text = MATH_D_RE.sub(lambda m: render_math(m.group(0)[2:-2], True, audit, name), text)
    text = MATH_I_RE.sub(lambda m: render_math(m.group(1), False, audit, name), text)
    audit["math_display"] += n_d
    audit["math_inline"] += n_i
    return text


# ---------------------------------------------------------------- code

# visibility policy (user, 2026-11-21): the docs show only the qmckl.h
# functions (C declarations: prototypes, typedefs, enums, defines) and ALL
# Fortran kernels; C function bodies, usage-example fragments, shell and
# plain-text blocks are hidden
C_BODY_RE = re.compile(r"\)\s*\{")  # a function definition / if/while body
COMPUTE_RE = re.compile(r"\bqmckl_compute_\w+\s*\(")
# the worked-example page shows complete programs, not API fragments
EXAMPLES_PAGE = "qmckl_examples.html"


def _c_nocomments(code: str) -> str:
    s = re.sub(r"/\*.*?\*/", "", code, flags=re.S)
    return re.sub(r"//[^\n]*", "", s)


def drop_reason(page: str, lang: str | None, code: str) -> str | None:
    """None = keep the block; otherwise why it is hidden (audit key)."""
    if page == EXAMPLES_PAGE:
        # complete programs: keep C, Fortran and shell blocks
        return None if lang in ("c", "C", "f90", "fortran", "bash") else "other-lang"
    if lang is None:
        return "plain"  # plain text / example-output blocks
    if lang in ("f90", "fortran"):
        return None
    if lang == "c":
        s = _c_nocomments(code)
        # qmckl_compute_* C prototypes are hidden too (user, 2026-11-21):
        # the substance of a compute function is its Fortran kernel
        if COMPUTE_RE.search(s):
            return "compute"
        if C_BODY_RE.search(s):
            return "body"
        return None
    return "other-lang"  # bash, emacs-lisp, makefile, ...


# ------------------------------------------------------------ worked examples
# Program sections of qmckl_examples.html -> downloadable file in code/.
# source "tangle": file is written by org-babel-tangle (export_docs.sh)
# source "blocks": verbatim concatenation of the section's f90 blocks
EXAMPLE_SECTIONS = [
    (r'<h3[^>]*>(?:<span[^>]*>[^<]*</span>)?\s*C</h3>',
     "mo_ortho.c", "tangle"),
    (r'<h3[^>]*>(?:<span[^>]*>[^<]*</span>)?\s*Fortran</h3>',
     "mo_ortho_f.f90", "blocks"),
    (r'<h2[^>]*>(?:<span[^>]*>[^<]*</span>)?\s*Checking errors</h2>',
     "qmckl_check_error.f90", "blocks"),
    (r'<h2[^>]*>(?:<span[^>]*>[^<]*</span>)?\s*'
     r'Computing an atomic orbital on a grid</h2>',
     "ao_grid.f90", "tangle"),
]


def example_programs(raw: str, audit: dict) -> tuple[dict, callable]:
    """Map every <pre> position of the examples page to its program file,
    and return {file: content} for the site's code/ directory."""
    sections = []  # (start, end, file)
    files = {}
    for hre, fname, source in EXAMPLE_SECTIONS:
        m = re.search(hre, raw)
        if m is None:
            audit["example_missing"].append(fname)
            continue
        dstart = raw.rfind('<div id="outline-container-', 0, m.start())
        _, _, dend = extract_balanced(raw, dstart)
        sections.append((dstart, dend, fname))
        if source == "blocks":
            blocks = [html.unescape(SPAN_RE.sub("", b)) for b in
                      re.findall(r'<pre class="src src-f90"[^>]*>(.*?)</pre>',
                                 raw[dstart:dend], re.S)]
            content = "\n".join(b for b in blocks if b.strip())
            if content and not content.endswith("\n"):
                content += "\n"
            files[fname] = content
    for hre, fname, source in EXAMPLE_SECTIONS:
        if source == "tangle":
            p = ORG / fname
            if not p.is_file() or not p.read_text().strip():
                raise SystemExit(
                    f"tangle file missing: {p} (run export_docs.sh first)")
            files[fname] = p.read_text()

    def file_at(pos: int) -> str | None:
        best = None
        for start, end, fname in sections:
            if start <= pos < end and (best is None or start > best[0]):
                best = (start, fname)
        return best[1] if best else None

    return files, file_at


def highlight_region(region: str, audit: dict, name: str,
                     file: str | None = None) -> str:
    """region = a full <pre ...>...</pre> (org source block)."""
    gt = region.index(">")
    tag, inner = region[: gt + 1], region[gt + 1 : -len("</pre>")]
    cls = re.search(r'class="([^"]*)"', tag)
    classes = cls.group(1).split() if cls else []
    lang = next((c[4:] for c in classes if c.startswith("src-")), None)
    audit["spans_stripped"] += len(SPAN_RE.findall(inner))
    code = html.unescape(SPAN_RE.sub("", inner))
    if not code.strip():
        # org emits an empty <pre> for src blocks whose content is
        # :export none (e.g. the "Provide" sections, generated by
        # ob-eval CALLs) — no empty card on the site
        audit["pre_empty"] += 1
        return ""
    reason = drop_reason(name, lang, code)
    if reason is not None:
        audit["pre_hidden"] += 1
        audit[f"pre_hidden_{reason}"] = audit.get(f"pre_hidden_{reason}", 0) + 1
        return ""
    label = ""
    if lang in LEXERS:
        lexer_name, label = LEXERS[lang]
        try:
            spans = highlight(code, get_lexer_by_name(lexer_name), HtmlFormatter(nowrap=True))
            audit["pre_highlighted"] += 1
        except Exception:
            spans = html.escape(code)
            audit["pre_fallback"].append({"file": name, "lang": lang})
    else:
        spans = html.escape(code)
        (audit["pre_unknown_lang"].append({"file": name, "lang": lang}) if lang
         else audit.__setitem__("pre_plain", audit["pre_plain"] + 1))
    # worked-example cards link to the complete, downloadable program file
    filelink = (f' <a class="code-file" href="../code/{file}" '
                f'title="Download {file}">{file}</a>') if file else ""
    cap = (f'<figcaption class="lang">{label}{filelink}</figcaption>'
           if (label or file) else "")
    fig = f'<figure class="code">{cap}' \
          f'<pre class="src-{lang or "plain"}"><code>{spans}</code></pre></figure>'
    return fig


def pass_code(text: str, regions: list[str], region_files: list[str | None],
              audit: dict, name: str) -> str:
    for i, region in enumerate(regions):
        replacement = (highlight_region(region, audit, name, region_files[i])
                       if region.startswith("<pre") else region)
        text = text.replace(f"\x00R{i:03d}\x00", replacement)
    return text


# ---------------------------------------------------------------- layout

def build_sidebar(current: str, toc_html: str) -> str:
    items = []
    for f in TOC:
        t = re.search(r"^#\+TITLE:\s*(.+)$", (ORG / f).read_text(), re.M)
        title = t.group(1).strip() if t else f.replace(".org", "").replace("_", " ")
        page = f.replace(".org", ".html")
        cur = ' aria-current="page"' if page == current else ""
        items.append(f'<li><a href="{page}"{cur}>{html.escape(title)}</a></li>')
    toc_part = f'  <nav class="side-toc" aria-label="On this page">\n{toc_html}\n  </nav>\n' if toc_html else ""
    return (
        '<aside class="sidebar">\n'
        '  <p class="side-brand"><a href="../index.html">QMCkl</a></p>\n'
        '  <nav class="side-nav" aria-label="Documentation">\n'
        f"    <ul>{''.join(items)}</ul>\n"
        "  </nav>\n"
        f"{toc_part}"
        "</aside>\n"
    )


def pass_layout(html_text: str, audit: dict, name: str) -> str:
    t = re.search(r"<title>(.*?)</title>", html_text, re.S)
    title = html.unescape(t.group(1)).strip() if t else name

    html_text = re.sub(r'<div id="org-div-home-and-up">.*?</div>', "", html_text, flags=re.S)
    html_text = re.sub(r'<div id="postamble".*?</div>\s*', "", html_text, flags=re.S)

    # page TOC -> sidebar
    toc_html = ""
    m = re.search(r'<div id="table-of-contents"[^>]*>', html_text)
    if m:
        c_start, content_end, region_end = extract_balanced(html_text, m.start())
        toc_block = html_text[c_start:content_end].replace("<h2>Table of Contents</h2>", "", 1)
        toc_html = f'<h2 class="side-toc-title">On this page</h2>{toc_block}'
        html_text = html_text[: m.start()] + html_text[region_end:]

    # page title: org emits either a title div or a bare <h1 class="title">
    h1 = ""
    mt = re.search(r'<div id="title"[^>]*>.*?</div>\s*|<h1 class="title"[^>]*>.*?</h1>\s*', html_text, re.S)
    if mt:
        mh = re.search(r"<h1[^>]*>(.*?)</h1>", mt.group(0), re.S)
        if mh:
            h1 = re.sub(r"\s+", " ", mh.group(1)).strip()
        html_text = html_text[: mt.start()] + html_text[mt.end():]

    # content div -> <main>
    mc = re.search(r'<div id="content"[^>]*>', html_text)
    if not mc:
        raise ValueError(f"{name}: no #content div")
    c_start, content_end, region_end = extract_balanced(html_text, mc.start())
    content = html_text[c_start:content_end]
    residual = (html_text[: mc.start()] + html_text[region_end:]).strip()
    if residual and residual != "":
        audit["residual_body"] = residual[:200]

    head = (
        "<head>\n"
        '<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f"<title>{html.escape(title)} — QMCkl</title>\n"
        f'<meta name="description" content="{html.escape(first_paragraph(content))}">\n'
        '<meta name="author" content="TREX CoE">\n'
        '<link rel="stylesheet" href="../css/theme.css">\n'
        "</head>"
    )
    return (
        '<!DOCTYPE html>\n<html lang="en">\n'
        f"{head}\n<body>\n"
        '<div class="layout">\n'
        f"{build_sidebar(name, toc_html)}\n"
        '<main id="content">\n'
        + (f'<header class="doc-head"><h1 class="title">{h1}</h1></header>\n' if h1 else "")
        + f"{content}\n"
        "</main>\n</div>\n</body>\n</html>\n"
    )


# ---------------------------------------------------------------- links

def fix_links(html_text: str, audit: dict, name: str) -> str:
    def sub(m):
        attr, quote, url = m.group(1), m.group(3), m.group(4)
        if url == "" or re.match(r"^(https?:|mailto:|#|tel:|data:)", url):
            return m.group(0)
        if url.endswith(".css"):
            return m.group(0)
        if url in ("index.html", "./index.html"):
            new = "../index.html"
        elif re.search(r"\.(png|jpe?g|svg|gif)$", url):
            new = f"../images/{url.split('/')[-1]}"
        else:
            new = url  # cross-doc link, same directory
        target = DOCS_OUT / new
        if not target.exists():
            audit["broken_links"].append({"file": name, "url": url})
        return f'{attr}={quote}{new}{quote}'

    return re.sub(r'(href|src)(=)(["\'])(.*?)\3', sub, html_text, flags=re.S)


# ---------------------------------------------------------------- table headers

def pass_row_headers(text: str, audit: dict, name: str) -> str:
    """Org tables without a header row (key/value and data-structure
    tables): make the first cell of each row a row header. Satisfies the
    axe td-has-header rule; the first column is the identity of the row
    (variable name, quantity, ...)."""
    out, pos = [], 0
    for tm in re.finditer(r'<table\b[^>]*>.*?</table>', text, re.S):
        seg = tm.group(0)
        if '<thead>' in seg:
            continue
        def rowconv(rm):
            row = rm.group(0)
            first = re.search(r'<td\b[^>]*>', row)
            if first:
                row = row[:first.start()] + '<th scope="row">' + row[first.end():]
                audit["row_headers"] += 1
            return row
        newseg = re.sub(r'<tr\b[^>]*>.*?</tr>', rowconv, seg, flags=re.S)
        out.append(text[pos:tm.start()])
        out.append(newseg)
        pos = tm.end()
    out.append(text[pos:])
    return "".join(out)


# ---------------------------------------------------------------- empty sections

def pass_empty_sections(text: str, audit: dict, name: str) -> str:
    """Drop outline sections left empty by :export none src blocks
    (e.g. the "Provide" sections: only ob-eval CALL results, no exported
    content) together with their sidebar TOC entries. Innermost levels
    first, so a parent is re-tested after its children are gone."""
    dropped = []
    for level in (4, 3, 2):
        opens = [m.start() for m in re.finditer(
            r'<div id="outline-container-[^"]*" class="outline-%d">' % level, text)]
        for start in reversed(opens):
            _, _, region_end = extract_balanced(text, start)
            inner = text[text.index(">", start) + 1: region_end - len("</div>")]
            # org HTML export: outline-N heading -> <hN> (outline-2 -> h2;
            # h1 is the page title, not an outline heading)
            hid = re.search(r'<h%d id="([^"]+)"' % level, inner)
            if hid is None:
                continue
            s = re.sub(r'<h%d\b.*?</h%d>' % (level, level), "", inner,
                       count=1, flags=re.S)
            s = re.sub(r'<div class="org-src-container">\s*</div>', "", s)
            txt = html.unescape(re.sub(r"<[^>]+>", " ", s))
            meaningful = ("<table", "<figure", "<math", "<img", "<ul", "<ol",
                          "<pre", "<blockquote", "<hr", 'math-display')
            # a section is empty only if its WHOLE remaining content (still
            # including children, which were re-tested at their own level in
            # an earlier pass) carries no visible content
            if (txt.strip() or any(t in s for t in meaningful)
                    or re.search(r'<div id="outline-container-', s)):
                continue
            text = text[:start] + text[region_end:]
            dropped.append(hid.group(1))
    for i in dropped:
        flat = r'<li><a href="#%s">.*?</a></li>' % re.escape(i)
        if re.search(flat, text):
            text = re.sub(flat, "", text, count=1)
            continue
        m = re.search(r'<li><a href="#%s">' % re.escape(i), text)
        if m:  # li with a nested <ul>: remove it balanced
            j, depth = m.start(), 0
            while True:
                mm = re.search(r"<li\b|</li>", text[j:])
                j += mm.end()
                depth += -1 if mm.group(0) == "</li>" else 1
                if depth == 0:
                    text = text[:m.start()] + text[j:]
                    break
    audit["empty_sections"] += len(dropped)
    return text


# ---------------------------------------------------------------- main

def process_doc(org_file: str, audit: dict) -> None:
    name = org_file.replace(".org", ".html")
    src = (ORG / name).read_text()
    a = {"file": name, "math_inline": 0, "math_display": 0, "math_eqnarray": 0,
         "math_fail": [],
         "pre_highlighted": 0, "pre_plain": 0, "pre_empty": 0, "pre_hidden": 0,
         "pre_fallback": [],
         "pre_unknown_lang": [], "spans_stripped": 0, "empty_sections": 0,
         "row_headers": 0,
         "example_files": {}, "example_missing": [],
         "broken_links": [], "residual_body": ""}

    # 0. worked examples: map each code region to its program file
    ex_files = {}
    file_at = (lambda pos: None)
    if name == EXAMPLES_PAGE:
        ex_files, file_at = example_programs(src, a)

    # 1. mask code regions (math and link fixing must not touch code)
    regions = []
    region_files = []
    def stash(m):
        regions.append(m.group(0))
        region_files.append(file_at(m.start()))
        return f"\x00R{len(regions) - 1:03d}\x00"
    text = PRE_RE.sub(stash, src)
    text = CODE_RE.sub(stash, text)

    # 2. math
    text = pass_math(text, a, name)
    # 3. layout (head, sidebar, main)
    text = pass_layout(text, a, name)
    # 4. links (code still masked)
    text = fix_links(text, a, name)
    # 5. restore code regions with highlighting
    text = pass_code(text, regions, region_files, a, name)
    # 6. tables without a header row: first cell -> row header (axe)
    text = pass_row_headers(text, a, name)
    # 7. drop sections left empty by hidden blocks (+ TOC entries)
    text = pass_empty_sections(text, a, name)

    (DOCS_OUT / name).write_text(text)
    if ex_files:
        code_dir = SITE / "code"
        code_dir.mkdir(exist_ok=True)
        for fname, content in ex_files.items():
            (code_dir / fname).write_text(content)
            a["example_files"][fname] = len(content.splitlines())
    audit[name] = a


def process_landing() -> None:
    """Highlight <pre data-lang="X"> blocks of tools/index.html and write
    OUT_DIR/index.html."""
    t = (TOOLSDIR / "index.html").read_text()

    def sub(m):
        lang, code = m.group(1), m.group(2)
        lexer_name, label = LEXERS.get(lang, (None, ""))
        try:
            spans = highlight(html.unescape(code),
                              get_lexer_by_name(lexer_name or "text"),
                              HtmlFormatter(nowrap=True))
        except Exception:
            spans = html.escape(code)
        return (
            f'<figure class="code">{f"<figcaption class=\"lang\">{label}</figcaption>" if label else ""}'
            f'<pre class="src-{lang}"><code>{spans}</code></pre></figure>'
        )

    t = re.sub(r'<pre data-lang="([^"]+)"[^>]*>(.*?)</pre>', sub, t, flags=re.S)
    (SITE / "index.html").write_text(t)


def main() -> None:
    DOCS_OUT.mkdir(parents=True, exist_ok=True)
    # landing page first: the docs link audit checks index.html on disk
    process_landing()
    audit = {}
    for f in TOC:
        process_doc(f, audit)

    report = {
        "pages": len(TOC),
        "math": {
            "inline": sum(a["math_inline"] for a in audit.values()),
            "display": sum(a["math_display"] for a in audit.values()),
            "eqnarray": sum(a["math_eqnarray"] for a in audit.values()),
            "failures": [x for a in audit.values() for x in a["math_fail"]],
        },
        "code": {
            "highlighted": sum(a["pre_highlighted"] for a in audit.values()),
            "plain": sum(a["pre_plain"] for a in audit.values()),
            "empty_blocks": sum(a["pre_empty"] for a in audit.values()),
            "hidden_blocks": sum(a["pre_hidden"] for a in audit.values()),
            "hidden_compute": sum(a.get("pre_hidden_compute", 0) for a in audit.values()),
            "spans_stripped": sum(a["spans_stripped"] for a in audit.values()),
            "fallback": [x for a in audit.values() for x in a["pre_fallback"]],
            "unknown_lang": [x for a in audit.values() for x in a["pre_unknown_lang"]],
        },
        "empty_sections": sum(a["empty_sections"] for a in audit.values()),
        "row_headers": sum(a["row_headers"] for a in audit.values()),
        "example_files": {f: n for a in audit.values() for f, n in a["example_files"].items()},
        "example_missing": [f for a in audit.values() for f in a["example_missing"]],
        "broken_links": [x for a in audit.values() for x in a["broken_links"]],
        "residual_body": {k: v for k, v in ((k, a["residual_body"]) for k, a in audit.items()) if v},
    }
    (SITE / "tools" / "_audit.json").write_text(json.dumps(report, indent=2))
    print(json.dumps({k: v for k, v in report.items() if k != "residual_body"}, indent=2))
    # executable <script> tags must be 0; <script type="application/ld+json"> is data, not JS
    scripts = sum(len(re.findall(r"<script(?![^>]*application/ld\+json)", p.read_text()))
                  for p in (SITE / "docs").glob("*.html"))
    scripts += len(re.findall(r"<script(?![^>]*application/ld\+json)", (SITE / "index.html").read_text()))
    print(f"script tags in output: {scripts} (must be 0)")


if __name__ == "__main__":
    main()
