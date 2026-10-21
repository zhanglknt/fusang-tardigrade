import zipfile, shutil, os, json, hashlib, csv

os.chdir(r"D:\系统发育树项目\Fusang\Fusang-main")
src = 'repro_package_v1.9.zip'
tmp = 'repro_v20_tmp'
assert not os.path.exists(tmp), 'tmp dir exists; clean it via shell first'
os.makedirs(tmp)
with zipfile.ZipFile(src) as z:
    z.extractall(tmp)

def sha256(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()

added = []

# 1. fastme_backend.py (imported by fusang_v2.py line ~934 for --tree_method fastme)
dst = os.path.join(tmp, 'fusang', 'fastme_backend.py')
shutil.copy2('fastme_backend.py', dst)
added.append(('fusang/fastme_backend.py', sha256('fastme_backend.py')))

# 2. generator script for the indel benchmark data
dst = os.path.join(tmp, 'scripts', 'gen_test_data_indel.py')
shutil.copy2('gen_test_data_indel.py', dst)
added.append(('scripts/gen_test_data_indel.py', sha256('gen_test_data_indel.py')))

# 3. nine per-seed CSVs referenced by DATA_SOURCES.json
csvs = ['benchmark_n200_clean_30seeds.csv', 'benchmark_n200_indel_30seeds.csv',
        'benchmark_n500_clean_30seeds.csv', 'benchmark_n500_indel_30seeds.csv',
        'benchmark_n1000_clean_30seeds.csv', 'table3_corrected.csv',
        'iqtree2_benchmark_indel_n200.csv', 'iqtree2_benchmark_indel_n200_full.csv',
        'benchmark_multik_ensemble_table7_full.csv']
for c in csvs:
    shutil.copy2(c, os.path.join(tmp, 'data', c))
    added.append((f'data/{c}', sha256(c)))

# 4. indel raw data seeds 100-129
raw = os.path.join(tmp, 'indel_raw_seeds100_129')
os.makedirs(raw)
for s in range(100, 130):
    for suffix in ('_indel.fasta', '_indel_true.nwk'):
        fn = f'seed{s}{suffix}'
        shutil.copy2(fn, os.path.join(raw, fn))
        added.append((f'indel_raw_seeds100_129/{fn}', sha256(fn)))
shutil.copy2('indel_raw_README.md', os.path.join(raw, 'README.md'))
added.append(('indel_raw_seeds100_129/README.md', sha256('indel_raw_README.md')))

# 5. DATA_SOURCES.json: record the additions
ds_path = os.path.join(tmp, 'data', 'DATA_SOURCES.json')
ds = json.load(open(ds_path, encoding='utf-8'))
ds['_v2.0_package'] = {
    'date': '2026-10-21',
    'note': ('Added fastme_backend.py (FastME wrapper imported by fusang_v2.py), '
             'all per-seed CSVs referenced by this manifest, the full 130-seed '
             'IQ-TREE2 GTR per-seed table (with per-seed status), and raw indel '
             'benchmark inputs for seeds 100-129 (FASTA + true trees) with the '
             'deterministic generator script gen_test_data_indel.py.'),
    'added_files': [{'path': p, 'sha256': h} for p, h in added],
}
json.dump(ds, open(ds_path, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)

# 6. README append
with open(os.path.join(tmp, 'README.md'), 'a', encoding='utf-8') as f:
    f.write('''

## v2.0 additions (2026-10-21)

- `fusang/fastme_backend.py` - FastME wrapper module imported by
  `fusang_v2.py` (`--tree_method fastme`); without it only the pure-Python
  NJ fallback runs
- `scripts/gen_test_data_indel.py` - the deterministic simulator that
  generated the indel benchmark data (JC69, sub=0.05, indel=0.02, L=500,
  geometric indel lengths mean=3; byte-identical regeneration given seed)
- `data/benchmark_{n200,n500,n1000}_{clean,indel}_30seeds.csv`,
  `data/table3_corrected.csv`, `data/iqtree2_benchmark_indel_n200*.csv`,
  `data/benchmark_multik_ensemble_table7_full.csv` - per-seed source data
  behind manuscript Tables 1/3/7 and the 130-seed benchmark (Figure 2)
- `indel_raw_seeds100_129/` - raw indel benchmark inputs: unaligned FASTA +
  true trees for seeds 100-129, with README documenting the generation path
''')

out = 'repro_package_v2.0.zip'
if os.path.exists(out):
    os.remove(out)
with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as zf:
    for root, dirs, files in os.walk(tmp):
        for fn in files:
            fp = os.path.join(root, fn)
            zf.write(fp, os.path.relpath(fp, tmp).replace(os.sep, '/'))
print('repro_package_v2.0.zip:', os.path.getsize(out), 'bytes,',
      len(zipfile.ZipFile(out).namelist()), 'files; added', len(added))
