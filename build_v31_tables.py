"""Rebuild All_Tables.docx for NAR submission v3.1 from NAR_MANUSCRIPT_REVISED.md (HEAD).

This is the v3.1 counterpart of step 3 in build_v30_package.py. It ONLY regenerates
All_Tables.docx -- it does not touch the manuscript DOCX, figures, or the repro ZIP.
"""
import os, re, shutil, subprocess, tempfile

SRC = r"D:\系统发育树项目\Fusang\Fusang-main"
MD = os.path.join(SRC, "NAR_MANUSCRIPT_REVISED.md")
PANDOC = r"C:\Users\admin\miniconda3\envs\human_selection\Library\bin\pandoc.exe"

with open(MD, encoding="utf-8") as f:
    content = f.read()

tables = re.findall(r'(\*\*Table (\d+)\. .+?\*\*)\n(.*?)(?=\n\n\*\*Table|\n\n##|\Z)',
                    content, re.DOTALL)
print(f"captured {len(tables)} tables: {[n for _, n, _ in tables]}")

tables_dir = os.path.join(SRC, "tables_pdf")
os.makedirs(tables_dir, exist_ok=True)
all_md = "\n\n\\newpage\n\n".join(f"{c.strip()}\n\n{b.strip()}" for c, n, b in tables)
tmp = os.path.join(tempfile.gettempdir(), "all_tables_v31.md")
open(tmp, "w", encoding="utf-8").write(all_md)

tables_docx = os.path.join(tables_dir, "All_Tables.docx")
r = subprocess.run([PANDOC, tmp, "-o", tables_docx, "--standalone"],
                   capture_output=True, text=True)
print(f"pandoc rc={r.returncode} {r.stderr[:200]}")
print(f"All_Tables.docx: {os.path.getsize(tables_docx):,} bytes")
os.remove(tmp)

# mirror into the v3.1 submission folder
dest_dir = os.path.join(SRC, "NAR_Submission_v3.1")
if os.path.isdir(dest_dir):
    shutil.copy2(tables_docx, os.path.join(dest_dir, "All_Tables.docx"))
    print(f"copied -> {dest_dir}\\All_Tables.docx")
