import zipfile, shutil, os

os.chdir(r"D:\系统发育树项目\Fusang\Fusang-main")
src = 'repro_package_v1.7.zip'
tmp = 'repro_v18_tmp'
if os.path.exists(tmp):
    shutil.rmtree(tmp)
os.makedirs(tmp)
with zipfile.ZipFile(src) as z:
    z.extractall(tmp)

GD = os.path.join('real_data', 'afproject_genome', 'gene_dna')
dstdir = os.path.join(tmp, 'data', 'afproject_genome', 'gene_dna')
os.makedirs(dstdir, exist_ok=True)
for item in os.listdir(GD):
    s = os.path.join(GD, item)
    t = os.path.join(dstdir, item)
    if os.path.isfile(s):
        shutil.copy2(s, t)
    elif os.path.isdir(s):
        os.makedirs(t, exist_ok=True)
        for fn in os.listdir(s):
            fp = os.path.join(s, fn)
            if os.path.isfile(fp) and os.path.getsize(fp) < 20_000_000:
                shutil.copy2(fp, os.path.join(t, fn))

# docs (v3.0 online documentation)
docs_dst = os.path.join(tmp, 'docs')
os.makedirs(docs_dst, exist_ok=True)
for f in ['TUTORIAL.md', 'EXAMPLES.md']:
    p = os.path.join('docs', f)
    if os.path.exists(p):
        shutil.copy2(p, os.path.join(docs_dst, f))
if os.path.exists('CITATION.cff'):
    shutil.copy2('CITATION.cff', os.path.join(tmp, 'CITATION.cff'))

with open(os.path.join(tmp, 'README.md'), 'a', encoding='utf-8') as f:
    f.write('''

## v1.8 additions (2026-10-09)

- `data/afproject_genome/gene_dna/` - gene-length real-DNA benchmark (manuscript
  v3.0, Supplementary Table S17): 13 mitochondrial protein-coding genes
  (168-1,866 bp) from the 25 fish genomes.
  - `gene_dna_results.json` - per-locus ETE3 nRF (vs Fischer et al. 2013
    reference and Fusang-vs-FastTree2 agreement), lengths, gap fractions
  - `Supplementary_Table_S17.md` - formatted supplementary table
  - `GENE_DNA_BENCHMARK.md` - methods + interpretation
  - `gbk/` (25 GenBank records), `fasta/` (13 multi-FASTA), `aln/` (MAFFT),
    `trees/` (per-locus Newick)
  - pipeline scripts: `fetch_gbk.py`, `extract_genes.py`, `run_fusang.py`,
    `run_mafft_ft2.sh`, `score_gene_dna.py`
  - key results: multi-k vs FT2 agreement nRF=0.451+-0.127; FT2 anchor vs
    reference 0.332+-0.175 (gene-tree discordance floor); multi-k 0.535 /
    k=7 0.532 / k=5 0.601 vs reference
- `docs/TUTORIAL.md`, `docs/EXAMPLES.md` - online documentation (v3.0);
  `CITATION.cff` - corrected authors and repository URL
''')

out = 'repro_package_v1.8.zip'
if os.path.exists(out):
    os.remove(out)
with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as zf:
    for root, dirs, files in os.walk(tmp):
        for fn in files:
            fp = os.path.join(root, fn)
            zf.write(fp, os.path.relpath(fp, tmp).replace(os.sep, '/'))
shutil.rmtree(tmp)
print('repro_package_v1.8.zip:', os.path.getsize(out), 'bytes,',
      len(zipfile.ZipFile(out).namelist()), 'files')
