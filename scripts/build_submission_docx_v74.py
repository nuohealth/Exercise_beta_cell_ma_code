#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build the DRCP submission DOCX for v4.0.

Single source of truth: manuscript/Manuscript_v6.3_DOM_20261029.md
Output: manuscript/Manuscript_v6.3_DOM_20261029.docx

Formatting applied (DRCP / Elsevier submission conventions):
  - A4, 2.5 cm margins
  - Times New Roman 12 pt, double line spacing, continuous line numbers
  - Tables as three-line (booktabs) tables, 9 pt, single spaced
  - Page numbers centred in footer
  - Figures appended at the end, one per page, with legends repeated
  - Internal build note ("> v3 change log") stripped from the journal-facing file

Run: python3 scripts/build_submission_docx_v3_1.py
"""
import os
import re
import subprocess
import sys

from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.enum.section import WD_SECTION
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(PROJ, "versions", "v7.4_20261031",
                   "Manuscript_v7.4_DOM_20261031.md")
BUILD_MD = os.path.join(PROJ, "manuscript", "_build_v74.md")
OUT = os.path.join(PROJ, "manuscript", "Manuscript_v7.4_DOM_20261031.docx")
FIGDIR = os.path.join(PROJ, "figures")

# 图片已按印刷宽度 16.0 cm 设计 → 插入宽度一律 16.0，1:1 不缩放
FIGS = [
    ("Fig1_PRISMA_v500.png", "Fig. 1", 16.0),
    ("Fig2_forest_three_axes_v500.png", "Fig. 2", 16.0),
    ("Fig3_decomposition_v600.png", "Fig. 3", 16.0),
    ("Fig4_DOSEEX_armwise.png", "Fig. 4", 16.0),
]


# ---------------------------------------------------------------- helpers
def strip_build_notes(text):
    """Remove internal-only blockquote lines (change logs, build notes)."""
    out = []
    for line in text.splitlines():
        if line.startswith("> "):
            continue
        out.append(line)
    return "\n".join(out)


def set_run_font(run, name="Times New Roman", size=12, bold=None):
    run.font.name = name
    run.font.size = Pt(size)
    if bold is not None:
        run.font.bold = bold
    rpr = run._element.get_or_add_rPr()
    rfonts = rpr.find(qn("w:rFonts"))
    if rfonts is None:
        rfonts = OxmlElement("w:rFonts")
        rpr.append(rfonts)
    for attr in ("w:ascii", "w:hAnsi", "w:cs", "w:eastAsia"):
        rfonts.set(qn(attr), name)


def add_line_numbering(section):
    """Continuous line numbering on every line."""
    sectPr = section._sectPr
    lnNumType = OxmlElement("w:lnNumType")
    lnNumType.set(qn("w:countBy"), "1")
    lnNumType.set(qn("w:restart"), "continuous")
    lnNumType.set(qn("w:distance"), "360")
    # schema order: ... pgMar, paperSrc, pgBorders, lnNumType, pgNumType, cols ...
    anchor = None
    for tag in ("w:pgNumType", "w:cols", "w:formProt", "w:vAlign"):
        el = sectPr.find(qn(tag))
        if el is not None:
            anchor = el
            break
    if anchor is not None:
        anchor.addprevious(lnNumType)
    else:
        sectPr.append(lnNumType)


def add_page_number_footer(section):
    """Single, well-formed PAGE field. Clears any pre-existing footer content first."""
    footer = section.footer
    footer.is_linked_to_previous = False
    p = footer.paragraphs[0] if footer.paragraphs else footer.add_paragraph()

    # wipe any existing runs / fields so the page number cannot be duplicated
    for child in list(p._p):
        if child.tag == qn("w:r"):
            p._p.remove(child)

    p.alignment = WD_ALIGN_PARAGRAPH.CENTER

    def _run():
        r = p.add_run()
        set_run_font(r, size=10)
        return r

    r1 = _run()
    fc_begin = OxmlElement("w:fldChar")
    fc_begin.set(qn("w:fldCharType"), "begin")
    r1._r.append(fc_begin)

    r2 = _run()
    instr = OxmlElement("w:instrText")
    instr.set(qn("xml:space"), "preserve")
    instr.text = " PAGE "
    r2._r.append(instr)

    r3 = _run()
    fc_sep = OxmlElement("w:fldChar")
    fc_sep.set(qn("w:fldCharType"), "separate")
    r3._r.append(fc_sep)

    r4 = _run()
    t = OxmlElement("w:t")
    t.text = "1"
    r4._r.append(t)

    r5 = _run()
    fc_end = OxmlElement("w:fldChar")
    fc_end.set(qn("w:fldCharType"), "end")
    r5._r.append(fc_end)


def three_line_table(table):
    """Booktabs-style: top rule, rule under header, bottom rule; no verticals."""
    tbl = table._tbl
    tblPr = tbl.tblPr
    for old in tblPr.findall(qn("w:tblBorders")):
        tblPr.remove(old)
    borders = OxmlElement("w:tblBorders")
    spec = {
        "top": ("single", "12"),
        "bottom": ("single", "12"),
        "left": ("none", "0"),
        "right": ("none", "0"),
        "insideH": ("none", "0"),
        "insideV": ("none", "0"),
    }
    for edge, (val, sz) in spec.items():
        el = OxmlElement(f"w:{edge}")
        el.set(qn("w:val"), val)
        if val != "none":
            el.set(qn("w:sz"), sz)
            el.set(qn("w:space"), "0")
            el.set(qn("w:color"), "000000")
        borders.append(el)
    tblPr.append(borders)

    # header row: rule beneath, bold, repeat across pages
    if len(table.rows):
        hdr = table.rows[0]
        trPr = hdr._tr.get_or_add_trPr()
        tblHeader = OxmlElement("w:tblHeader")
        tblHeader.set(qn("w:val"), "true")
        trPr.append(tblHeader)
        for cell in hdr.cells:
            tcPr = cell._tc.get_or_add_tcPr()
            tcBorders = OxmlElement("w:tcBorders")
            bottom = OxmlElement("w:bottom")
            bottom.set(qn("w:val"), "single")
            bottom.set(qn("w:sz"), "8")
            bottom.set(qn("w:space"), "0")
            bottom.set(qn("w:color"), "000000")
            tcBorders.append(bottom)
            tcPr.append(tcBorders)

    # single spacing + 9 pt everywhere in the table
    for row in table.rows:
        for cell in row.cells:
            for p in cell.paragraphs:
                pf = p.paragraph_format
                pf.line_spacing_rule = WD_LINE_SPACING.SINGLE
                pf.space_before = Pt(1)
                pf.space_after = Pt(1)
                for r in p.runs:
                    set_run_font(r, size=9)


def set_col_widths(table, widths_cm, font_pt=None):
    """Give each column an explicit width so headers do not break mid-word."""
    tbl = table._tbl
    tblPr = tbl.tblPr

    # layout = fixed so the widths we set are honoured
    for old in tblPr.findall(qn("w:tblLayout")):
        tblPr.remove(old)
    layout = OxmlElement("w:tblLayout")
    layout.set(qn("w:type"), "fixed")
    tblPr.append(layout)

    # total table width in twentieths of a point
    total_tw = int(sum(widths_cm) * 567)
    for old in tblPr.findall(qn("w:tblW")):
        tblPr.remove(old)
    tblW = OxmlElement("w:tblW")
    tblW.set(qn("w:w"), str(total_tw))
    tblW.set(qn("w:type"), "dxa")
    tblPr.append(tblW)

    # grid columns
    grid = tbl.find(qn("w:tblGrid"))
    if grid is not None:
        for gc, w in zip(grid.findall(qn("w:gridCol")), widths_cm):
            gc.set(qn("w:w"), str(int(w * 567)))

    for row in table.rows:
        for cell, w in zip(row.cells, widths_cm):
            cell.width = Cm(w)
            tcPr = cell._tc.get_or_add_tcPr()
            for old in tcPr.findall(qn("w:tcW")):
                tcPr.remove(old)
            tcW = OxmlElement("w:tcW")
            tcW.set(qn("w:w"), str(int(w * 567)))
            tcW.set(qn("w:type"), "dxa")
            tcPr.append(tcW)

    if font_pt:
        for row in table.rows:
            for cell in row.cells:
                for para in cell.paragraphs:
                    for r in para.runs:
                        set_run_font(r, size=font_pt)


# Column widths (cm). Total must stay <= 16.0 (A4 minus 2.5 cm margins).
TABLE_WIDTHS = {
    # Table 1, 10 cols: Reference Study Design n Duration Baseline glucose Exercise Secretion measure Provenance Clearance
    # 注：表在 v3.13 加了 Baseline glucose 列后变 10 列，旧表只有 9 列宽 → len 不匹配，
    # set_col_widths 被整段跳过，全部列退化成均宽，表头断词（Referen/Duratio）。v3.15 修正。
    10: ([1.22, 1.98, 1.40, 0.70, 1.28, 2.42, 2.05, 2.05, 1.55, 1.35], 6.8),
    # Table 2, 9 cols: Axis Stratum k(arms) k(studies) Pooled change 95% CI p value I2 tau2
    # 注：v3.13 起 Table 2 由 7 列扩为 9 列，旧的 7 列宽失效 → 退化成均宽，表头断词（Dispositi/on、sensitivit/y）。v3.15 修正。
    9: ([2.80, 2.05, 1.25, 1.35, 1.85, 2.50, 1.50, 1.15, 1.10], 8.5),
    # 旧版 7 列（保留兼容，若某表回落到 7 列仍可用）
    7: ([3.30, 1.95, 1.95, 1.95, 1.95, 1.95, 1.95], 9.0),
    # Table 3, 8 cols: Axis k(arms) k(studies) Pooled 95%CI p I2 tau2
    8: ([3.60, 1.45, 1.60, 2.10, 3.10, 1.40, 1.35, 1.40], 8.5),
}


def tighten_cell_margins(table, left_cm=0.06, right_cm=0.06):
    """Reduce default cell padding so wide tables fit the text block."""
    tblPr = table._tbl.tblPr
    for old in tblPr.findall(qn("w:tblCellMar")):
        tblPr.remove(old)
    mar = OxmlElement("w:tblCellMar")
    for side, cm in (("left", left_cm), ("right", right_cm)):
        el = OxmlElement(f"w:{side}")
        el.set(qn("w:w"), str(int(cm * 567)))
        el.set(qn("w:type"), "dxa")
        mar.append(el)
    tblPr.append(mar)


def keep_table_together(table):
    """Stop a table from stranding its last rows on the next page.

    A table this short should never break: the row-level "do not split" flag
    prevents a single row being halved, and keep-with-next on every row but the
    last keeps the whole block on one page. Without this, Table 2 in v3.20 put
    its last two rows on the page that opens with the Table 3 title.
    """
    rows = table.rows
    for i, row in enumerate(rows):
        trPr = row._tr.get_or_add_trPr()
        if trPr.find(qn("w:cantSplit")) is None:
            trPr.append(OxmlElement("w:cantSplit"))
        if i < len(rows) - 1:
            for cell in row.cells:
                for p in cell.paragraphs:
                    p.paragraph_format.keep_with_next = True


def apply_table_layout(doc):
    for table in doc.tables:
        ncol = len(table.columns)
        spec = TABLE_WIDTHS.get(ncol)
        if spec:
            widths, pt = spec
            # ensure the width list matches the actual column count
            if len(widths) == ncol:
                set_col_widths(table, widths, pt)
                tighten_cell_margins(table)


def style_body(doc):
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(12)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    pf = st.paragraph_format
    pf.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    pf.space_before = Pt(0)
    pf.space_after = Pt(0)

    for name, size in (("Heading 1", 14), ("Heading 2", 13), ("Heading 3", 12), ("Title", 16)):
        try:
            s = doc.styles[name]
        except KeyError:
            continue
        s.font.name = "Times New Roman"
        s.font.size = Pt(size)
        s.font.bold = True
        s.font.color.rgb = RGBColor(0, 0, 0)
        s.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
        s.paragraph_format.space_before = Pt(12)
        s.paragraph_format.space_after = Pt(6)
        s.paragraph_format.keep_with_next = True


# ---------------------------------------------------------------- superscripts
CITE_RE = re.compile(r'\[(\d+(?:[,\u2013\-]\d+)*)\]')


def superscript_citations(doc):
    """DOM: in-text citations must be Arabic superscript numerals, not [n]."""
    # find the References heading; everything after it is left alone
    stop = len(doc.paragraphs)
    for i, para in enumerate(doc.paragraphs):
        if para.text.strip() == "References":
            stop = i
            break

    def convert(para):
        if not any(CITE_RE.search(r.text or "") for r in para.runs):
            return 0
        # merge runs' text, then rebuild preserving per-run formatting
        new_runs = []
        n = 0
        for run in list(para.runs):
            txt = run.text or ""
            if not CITE_RE.search(txt):
                new_runs.append((txt, False, run))
                continue
            pos = 0
            for m in CITE_RE.finditer(txt):
                if m.start() > pos:
                    new_runs.append((txt[pos:m.start()], False, run))
                new_runs.append((m.group(1), True, run))
                n += 1
                pos = m.end()
            if pos < len(txt):
                new_runs.append((txt[pos:], False, run))
        if n == 0:
            return 0
        # wipe and re-add
        for run in list(para.runs):
            run._r.getparent().remove(run._r)
        for txt, sup, src in new_runs:
            if txt == "":
                continue
            r = para.add_run(txt)
            r.bold = src.bold
            r.italic = src.italic
            if sup:
                r.font.superscript = True
        return n

    total = 0
    for para in doc.paragraphs[:stop]:
        total += convert(para)
    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                for para in cell.paragraphs:
                    total += convert(para)
    return total


def main():
    if not os.path.exists(SRC):
        sys.exit(f"source not found: {SRC}")

    text = open(SRC, encoding="utf-8").read()
    text = strip_build_notes(text)
    # journal-facing header line: normalise to the current version stamp
    text = re.sub(
        r"^Manuscript v[0-9.]+ — DRCP format — \d{4}-\d{2}-\d{2}$",
        "Manuscript v7.0 — DOM format — 2026-10-31",
        text, count=1, flags=re.MULTILINE,
    )
    open(BUILD_MD, "w", encoding="utf-8").write(text)

    subprocess.run(
        ["pandoc", BUILD_MD, "-f", "markdown+pipe_tables-yaml_metadata_block", "-o", OUT, "--columns=1000"],
        check=True,
    )

    doc = Document(OUT)
    style_body(doc)

    for section in doc.sections:
        section.page_width = Cm(21.0)
        section.page_height = Cm(29.7)
        section.top_margin = Cm(2.5)
        section.bottom_margin = Cm(2.5)
        section.left_margin = Cm(2.5)
        section.right_margin = Cm(2.5)
        add_line_numbering(section)
        add_page_number_footer(section)

    for table in doc.tables:
        three_line_table(table)

    n_sup = superscript_citations(doc)
    print(f"  superscript citations: {n_sup}")

    apply_table_layout(doc)

    for table in doc.tables:
        keep_table_together(table)

    # figures at the very end, one per page, legend repeated
    doc.add_page_break()
    p = doc.add_paragraph()
    r = p.add_run("Figures")
    set_run_font(r, size=14, bold=True)

    for fname, label, width_cm in FIGS:
        path = os.path.join(FIGDIR, fname)
        if not os.path.exists(path):
            print(f"  ! missing figure: {path}")
            continue
        doc.add_page_break()
        cap = doc.add_paragraph()
        rr = cap.add_run(label)
        set_run_font(rr, size=11, bold=True)
        doc.add_picture(path, width=Cm(width_cm))
        doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER

    doc.save(OUT)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
