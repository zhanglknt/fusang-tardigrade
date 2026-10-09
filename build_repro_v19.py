import zipfile, shutil, os

os.chdir(r"D:\系统发育树项目\Fusang\Fusang-main")
src = 'repro_package_v1.8.zip'
tmp = 'repro_v19_tmp'
if os.path.exists(tmp):
    shutil.rmtree(tmp)
os.makedirs(tmp)
with zipfile.ZipFile(src) as z:
    z.extractall(tmp)

# FastME binaries (default tree builder for simplified pipeline)
fm_dst = os.path.join(tmp, 'fusang', 'fastme_bin')
os.makedirs(fm_dst, exist_ok=True)
for f in ['fastme.exe', 'fastme_linux']:
    shutil.copy2(os.path.join('fastme_bin', f), os.path.join(fm_dst, f))

with open(os.path.join(tmp, 'README.md'), 'a', encoding='utf-8') as f:
    f.write('''

## v1.9 additions (2026-10-09)

- `fusang/fastme_bin/fastme.exe`, `fusang/fastme_bin/fastme_linux` - FastME
  v2.1.6.4 binaries (Windows native + Linux), the default tree builder
  (BIONJ+BNNI) for the simplified pipeline; fusang_v2.py discovers them via
  its priority cascade (bundled binary -> system fastme -> BioPython NJ
  fallback)
''')

out = 'repro_package_v1.9.zip'
if os.path.exists(out):
    os.remove(out)
with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as zf:
    for root, dirs, files in os.walk(tmp):
        for fn in files:
            fp = os.path.join(root, fn)
            zf.write(fp, os.path.relpath(fp, tmp).replace(os.sep, '/'))
shutil.rmtree(tmp)
print('repro_package_v1.9.zip:', os.path.getsize(out), 'bytes,',
      len(zipfile.ZipFile(out).namelist()), 'files')
