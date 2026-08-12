#!/usr/bin/env python
"""Render plans/active/titan-campaign-methods.md -> .html.

Docs ruling (task #103): the .html must NOT be a hand-maintained duplicate of
the .md (drift risk). This makes the .md the single source of truth and
regenerates the .html from it, stamping a generator marker so no one edits the
.html by hand. Re-run after any edit to the .md:

    mamba run -n PPcl python plans/scripts/render_methods_html.py

The bespoke print/paper CSS shell is preserved here; the <body> is generated
from the markdown (tables, fenced code, and the leading '---' rule supported).
"""
from pathlib import Path

import markdown

REPO = Path(__file__).resolve().parents[2]
MD = REPO / 'plans/active/titan-campaign-methods.md'
HTML = REPO / 'plans/active/titan-campaign-methods.html'

STYLE = """  :root {
    --fg: #1a1a1a; --muted: #555; --rule: #d8d8d8;
    --accent: #1f5c8f; --band: #f4f6f8; --code: #f0f1f3;
    --warn-bg: #fdf3e7; --warn-br: #d98324; --ok: #2b7a4b;
  }
  html { -webkit-text-size-adjust: 100%; }
  body {
    font-family: "Iowan Old Style", "Palatino Linotype", Palatino, Georgia, serif;
    color: var(--fg); line-height: 1.55; max-width: 820px;
    margin: 0 auto; padding: 3rem 1.5rem 6rem; font-size: 17px;
  }
  h1 { font-size: 1.9rem; line-height: 1.2; margin: 0 0 .3rem; }
  h2 {
    font-size: 1.3rem; margin: 2.6rem 0 .8rem; padding-bottom: .25rem;
    border-bottom: 2px solid var(--accent); color: var(--accent);
  }
  h3 { font-size: 1.05rem; margin: 1.6rem 0 .4rem; }
  p, li { margin: .5rem 0; }
  code, .mono {
    font-family: "SF Mono", "JetBrains Mono", Menlo, Consolas, monospace;
    font-size: .86em; background: var(--code); padding: .05em .35em;
    border-radius: 3px;
  }
  pre {
    background: var(--band); border: 1px solid var(--rule); border-radius: 6px;
    padding: .9rem 1.1rem; overflow-x: auto; font-size: .85rem; line-height: 1.5;
  }
  pre code { background: none; padding: 0; }
  table {
    border-collapse: collapse; width: 100%; margin: 1rem 0; font-size: .9rem;
  }
  th, td { text-align: left; padding: .4rem .6rem; border-bottom: 1px solid var(--rule); }
  th { border-bottom: 2px solid var(--accent); color: var(--accent); font-weight: 600; }
  tbody tr:nth-child(odd) { background: var(--band); }
  strong { color: var(--fg); }
  hr { border: none; border-top: 1px solid var(--rule); margin: 2rem 0; }
  sub { font-size: .75em; }
  blockquote {
    border-left: 4px solid var(--warn-br); background: var(--warn-bg);
    padding: .8rem 1.1rem; margin: 1.2rem 0; border-radius: 0 6px 6px 0;
  }"""

MARKER = ('<!-- GENERATED FILE — DO NOT EDIT BY HAND.\n'
          '     Source: plans/active/titan-campaign-methods.md\n'
          '     Regenerate: mamba run -n PPcl python '
          'plans/scripts/render_methods_html.py -->')


def main():
    body = markdown.markdown(
        MD.read_text(),
        extensions=['tables', 'fenced_code', 'sane_lists', 'attr_list'],
    )
    doc = (
        '<!DOCTYPE html>\n<html lang="en">\n<head>\n'
        '<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f'{MARKER}\n'
        '<title>Titan free-gravity SBI campaign — methods</title>\n'
        f'<style>\n{STYLE}\n</style>\n</head>\n<body>\n\n'
        f'{body}\n\n</body>\n</html>\n'
    )
    HTML.write_text(doc)
    print(f'wrote {HTML} ({len(doc)} chars) from {MD.name}')


if __name__ == '__main__':
    main()
