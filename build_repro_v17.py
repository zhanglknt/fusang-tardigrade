import zipfile, shutil, os

os.chdir(r"D:\系统发育树项目\Fusang\Fusang-main")
src = 'repro_package_v1.5.zip'
tmp = 'repro_v16_tmp'
if os.path.exists(tmp):
    shutil.rmtree(tmp)
os.makedirs(tmp)
with zipfile.ZipFile(src) as z:
    z.extractall(tmp)

AG = os.path.join('real_data', 'afproject_genome')
dstdir = os.path.join(tmp, 'data', 'afproject_genome')
os.makedirs(dstdir, exist_ok=True)
for f in ['real_genome_scores.json', 'Supplementary_Table_S16.md',
          'fish_mito_suite.py', 'ecoli_shigella_suite.py', 'ete3_scoring.py',
          'score_external.py', 'run_real_genome_wsl.sh', 'swisstree_gap_stats.sh']:
    shutil.copy2(os.path.join(AG, f), os.path.join(dstdir, f))

fishdir = os.path.join(dstdir, 'fish_mito')
os.makedirs(fishdir, exist_ok=True)
for fn in os.listdir(os.path.join(AG, 'fish_mito')):
    if fn.endswith('.fasta'):
        shutil.copy2(os.path.join(AG, 'fish_mito', fn), os.path.join(fishdir, fn))

for sub in ['fish_mito_results', 'ecoli_results', 'ecoli_results_canon', 'wsl_out']:
    s = os.path.join(AG, sub)
    t = os.path.join(dstdir, sub)
    if os.path.isdir(s):
        os.makedirs(t, exist_ok=True)
        for fn in os.listdir(s):
            fp = os.path.join(s, fn)
            if os.path.isfile(fp) and os.path.getsize(fp) < 20_000_000:
                shutil.copy2(fp, os.path.join(t, fn))

with open(os.path.join(tmp, 'README.md'), 'a', encoding='utf-8') as f:
    f.write('''

## v1.7 additions (2026-10-09)

- `data/afproject_genome/` - real indel-rich genome benchmarks (AFproject community
  datasets, reservation 3 of the review closure):
  - `real_genome_scores.json` - unified ETE3 nRF results, all methods, both datasets
  - `Supplementary_Table_S16.md` - full supplementary table (accessions, gap stats,
    per-method results, published-value cross-validation)
  - `fish_mito/` - 25 fish mitochondrial genomes (NCBI); `fish_mito_results/` - trees
  - `ecoli_results/` / `ecoli_results_canon/` - E. coli/Shigella genome-cosine trees
    (forward-strand and canonical variants)
  - `wsl_out/` - kmacs dmat / Mash tsv distance matrices, MAFFT alignment, FastTree2
    tree; E. coli genome FASTAs are NOT bundled (130 MB) - download via NCBI
    accessions listed in Supplementary_Table_S16.md
  - suite scripts: `fish_mito_suite.py`, `ecoli_shigella_suite.py` (64-bit rolling
    hash k-mer cosine), `ete3_scoring.py`, `run_real_genome_wsl.sh`,
    `swisstree_gap_stats.sh`
  - key results: fish multi-k cosine nRF=0.045 (ties MSA+ML and best published AF
    methods); E. coli/Shigella boundary (cosine 0.50 vs Mash 0.21/0.12 published,
    anchors 0.08 published); SwissTree gap content median 50%
''')

out = 'repro_package_v1.7.zip'
if os.path.exists(out):
    os.remove(out)
with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as zf:
    for root, dirs, files in os.walk(tmp):
        for fn in files:
            fp = os.path.join(root, fn)
            zf.write(fp, os.path.relpath(fp, tmp).replace(os.sep, '/'))
shutil.rmtree(tmp)
print('repro_package_v1.7.zip:', os.path.getsize(out), 'bytes,',
      len(zipfile.ZipFile(out).namelist()), 'files')
