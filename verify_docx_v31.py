"""Self-check for NAR_MANUSCRIPT_REVISED_v3.1.docx.

Reopens the generated docx, reports structure statistics, and greps for
current (must be present) and superseded (must be judged case by case) values.
"""
import os
import re
import sys

from docx import Document

SRC = os.path.dirname(os.path.abspath(__file__))
DOCX = os.path.join(SRC, "NAR_MANUSCRIPT_REVISED_v3.1.docx")
MD = os.path.join(SRC, "NAR_MANUSCRIPT_REVISED.md")

MUST_PRESENT = ["0.0507", "2.07\u00d710\u207b\u00b9\u2074", "0.669", "0.2617",
                "0.2726", "1.54\u00d710\u207b\u2074", "0.294"]
WATCH = ["p=0.049", "2.3\u00d710\u207b\u00b9\u2074", "0.709", "0.083 \u00b1 0.005",
         "regardless of scale", "indel rates above 0.01"]

# p-value metric switch (float nRF -> integer RF), applied by stats-align
P_NEW = ["2.22\u00d710\u207b\u00b9\u00b2", "9.21\u00d710\u207b\u2075",
         "0.00206", "0.00247"]
P_OLD = ["2.68\u00d710\u207b\u00b9\u00b2", "7.95\u00d710\u207b\u2075",
         "0.0018", "0.0024"]

SECTIONS = ["ABSTRACT", "INTRODUCTION", "MATERIALS AND METHODS", "RESULTS",
            "DISCUSSION", "SUPPLEMENTARY MATERIAL", "DATA AVAILABILITY",
            "FUNDING", "COMPETING INTERESTS", "ACKNOWLEDGEMENTS",
            "AUTHOR CONTRIBUTIONS", "CORRESPONDING AUTHOR", "FIGURE LEGENDS",
            "REFERENCES"]


def count_exact(text, value):
    """Count occurrences that are not part of a longer number
    (so that 0.0024 does not match inside 0.00247)."""
    return len(re.findall(re.escape(value) + r"(?![0-9])", text))


def table_blocks(doc):
    """For each table: its cell text plus the footnote paragraphs that follow it
    (up to the next table or the next section heading)."""
    from docx.table import Table
    from docx.text.paragraph import Paragraph

    items = []
    for child in doc.element.body.iterchildren():
        tag = child.tag.split("}")[-1]
        if tag == "p":
            items.append(("p", Paragraph(child, doc)))
        elif tag == "tbl":
            items.append(("t", Table(child, doc)))

    blocks = []
    for k, (kind, obj) in enumerate(items):
        if kind != "t":
            continue
        parts = ["\n".join(c.text for r in obj.rows for c in r.cells)]
        for kind2, obj2 in items[k + 1:]:
            if kind2 == "t":
                break
            txt = obj2.text.strip()
            if txt in SECTIONS:
                break
            if txt:
                parts.append(txt)
            if len(parts) > 6:
                break
        blocks.append("\n".join(parts))
    return blocks


def section_text(doc, name):
    """Text of a top-level section (heading paragraph up to the next heading)."""
    idx = [i for i, p in enumerate(doc.paragraphs) if p.text.strip() in SECTIONS]
    pos = [i for i in idx if doc.paragraphs[i].text.strip() == name]
    if not pos:
        return ""
    start = pos[0]
    end = next((i for i in idx if i > start), len(doc.paragraphs))
    return "\n".join(p.text for p in doc.paragraphs[start:end])


def full_text(doc):
    parts = []
    for p in doc.paragraphs:
        parts.append(p.text)
    for t in doc.tables:
        for row in t.rows:
            for cell in row.cells:
                parts.append(cell.text)
    return "\n".join(parts)


def main():
    doc = Document(DOCX)
    text = full_text(doc)

    n_par = len(doc.paragraphs)
    n_tab = len(doc.tables)
    rows = [len(t.rows) for t in doc.tables]
    n_rows = sum(rows)
    n_sec = len(doc.sections)

    print("=" * 68)
    print("STRUCTURE")
    print("=" * 68)
    print(f"paragraphs (top level)   : {n_par}")
    print(f"tables                   : {n_tab}   (expected 12)")
    print(f"table rows               : {n_rows}")
    print(f"per-table rows           : {rows}")
    print(f"per-table cols           : {[len(t.columns) for t in doc.tables]}")
    print(f"sections                 : {n_sec}")
    print(f"page size                : {doc.sections[0].page_width.cm:.1f} x "
          f"{doc.sections[0].page_height.cm:.1f} cm")

    words = len(re.findall(r"\S+", text))
    size_kb = os.path.getsize(DOCX) / 1024
    print(f"word count               : {words}")
    print(f"file size                : {size_kb:.1f} KB")
    print(f"estimated pages (w/600)  : {words / 600:.1f}")

    print()
    print("=" * 68)
    print("MUST BE PRESENT")
    print("=" * 68)
    for v in MUST_PRESENT:
        c = text.count(v)
        print(f"  [{'OK ' if c else 'MISS'}] {v!r}: {c} occurrence(s)")

    print()
    print("=" * 68)
    print("SUPERSEDED VALUES - CONTEXT INSPECTION")
    print("=" * 68)
    for v in WATCH:
        hits = [m.start() for m in re.finditer(re.escape(v), text)]
        print(f"\n  {v!r}: {len(hits)} hit(s)")
        for h in hits:
            ctx = text[max(0, h - 170):h + 170].replace("\n", " ")
            print(f"     ...{ctx}...")

    print()
    print("=" * 68)
    print("MD <-> DOCX COVERAGE (every non-table source line must appear)")
    print("=" * 68)
    with open(MD, encoding="utf-8") as fh:
        md_lines = fh.read().split("\n")

    def flat(s):
        # stash `code` spans so that '*' inside them is not stripped
        spans = []

        def _stash(m):
            spans.append(m.group(1))
            return "\x00%d\x00" % (len(spans) - 1)

        s = re.sub(r"`([^`]+)`", _stash, s)
        s = re.sub(r"[*`#]", "", s)
        s = re.sub(r"\x00(\d+)\x00", lambda m: spans[int(m.group(1))], s)
        s = re.sub(r"^[\-\u2022]\s+", "", s.strip())
        return re.sub(r"\s+", " ", s).strip()

    flat_docx = re.sub(r"\s+", " ", text)
    md_para = [l.strip() for l in md_lines
               if l.strip() and l.strip() != "---" and not l.strip().startswith("|")]
    missing = [l for l in md_para if flat(l) and flat(l) not in flat_docx]
    print(f"md non-table lines        : {len(md_para)}")
    print(f"non-table lines not found : {len(missing)}")
    for l in missing[:20]:
        print("   MISSING:", l[:160])

    # markdown table cell coverage
    md_cells = []
    for l in md_lines:
        if l.strip().startswith("|") and not re.match(r"^\|[\s:|-]+\|$", l.strip()):
            md_cells += [c.strip() for c in l.strip().strip("|").split("|")]
    cells_missing = [c for c in md_cells if c and flat(c) not in flat_docx]
    print(f"md table cells            : {len(md_cells)}")
    print(f"table cells not found     : {len(cells_missing)}")
    for c in cells_missing[:20]:
        print("   MISSING CELL:", c[:120])

    print()
    print("=" * 68)
    print("REFERENCE NUMBERING (REFERENCES section only)")
    print("=" * 68)
    m = re.search(r"REFERENCES(.*?)(?:FIGURE LEGENDS|Manuscript prepared)", text, re.S)
    seg = m.group(1) if m else ""
    refs = re.findall(r"^\s*(\d+)\.\s+[A-Z]", seg, re.MULTILINE)
    print("reference entries         :", len(refs))
    print("numbering 1..N in order   :",
          refs == [str(k) for k in range(1, len(refs) + 1)])
    print("numbers                   :", refs)

    print()
    print("=" * 68)
    print("P-VALUE SWITCH (integer-RF restatement by stats-align)")
    print("=" * 68)
    for v in P_NEW:
        print(f"  NEW {v!r}: {count_exact(text, v)} occurrence(s)")
    for v in P_OLD:
        hits = [m.start() for m in
                re.finditer(re.escape(v) + r"(?![0-9])", text)]
        print(f"\n  OLD {v!r}: {len(hits)} hit(s)")
        for h in hits:
            print("     ..." + text[max(0, h - 150):h + 150].replace("\n", " ") + "...")

    print()
    print("=" * 68)
    print("TARGETED CHECK 1: Table 4 / Table 5 cells + their footnotes")
    print("=" * 68)
    blocks = table_blocks(doc)  # cells + paragraphs trailing each table
    for ti in (3, 4):  # 0-based: Table 4 -> idx 3, Table 5 -> idx 4
        if ti >= len(blocks):
            continue
        blob = blocks[ti]
        print(f"  table index {ti} (Table {ti + 1}):")
        found = False
        for v in P_NEW + P_OLD:
            c = count_exact(blob, v)
            if c:
                found = True
                print(f"     {v!r}: {c}")
                h = re.search(re.escape(v) + r"(?![0-9])", blob).start()
                print("        ..." + blob[max(0, h - 120):h + 120].replace("\n", " ") + "...")
        if not found:
            print("     (no p-value switch values present)")

    print()
    print("=" * 68)
    print("TARGETED CHECK 2: ABSTRACT section")
    print("=" * 68)
    abst = section_text(doc, "ABSTRACT")
    print(f"  ABSTRACT length: {len(abst)} chars")
    for v in P_NEW + P_OLD:
        c = count_exact(abst, v)
        if c:
            h = re.search(re.escape(v) + r"(?![0-9])", abst).start()
            print(f"  {v!r}: {c} -> ..." + abst[max(0, h - 110):h + 110] + "...")

    print()
    print("=" * 68)
    print("SECTION ORDER IN DOCX")
    print("=" * 68)
    order = [p.text.strip() for p in doc.paragraphs
             if re.match(r"^(ABSTRACT|INTRODUCTION|MATERIALS AND METHODS|RESULTS|"
                         r"DISCUSSION|SUPPLEMENTARY MATERIAL|DATA AVAILABILITY|"
                         r"FUNDING|COMPETING INTERESTS|ACKNOWLEDGEMENTS|"
                         r"AUTHOR CONTRIBUTIONS|CORRESPONDING AUTHOR|REFERENCES|"
                         r"FIGURE LEGENDS)$", p.text.strip())]
    print(" -> ".join(order))


if __name__ == "__main__":
    main()
    sys.exit(0)
