#!/usr/bin/env python3
"""Step 1: fetch GenBank records for 25 fish_mito accessions via NCBI eutils."""
import os, sys, time, urllib.request

BASE = os.path.dirname(os.path.abspath(__file__))
GBK = os.path.join(BASE, "gbk")
os.makedirs(GBK, exist_ok=True)

DATA = r"D:\系统发育树项目\Fusang\Fusang-main\real_data\afproject_genome\fish_mito"
accs = [fn[:-6] for fn in sorted(os.listdir(DATA)) if fn.endswith(".fasta")]
print(f"{len(accs)} accessions")

UA = "fusang-bench/1.0 (mailto:knightz@pumc.edu.cn)"
for acc in accs:
    out = os.path.join(GBK, acc + ".gbk")
    if os.path.exists(out) and os.path.getsize(out) > 1000:
        print(f"  {acc}: cached"); continue
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={acc}&rettype=gb&retmode=text")
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = r.read()
            if b"LOCUS" not in data:
                raise RuntimeError("not a genbank record")
            open(out, "wb").write(data)
            print(f"  {acc}: {len(data)} bytes")
            break
        except Exception as e:
            print(f"  {acc}: attempt {attempt+1} failed: {e}", file=sys.stderr)
            time.sleep(2)
    else:
        print(f"  {acc}: FAILED"); continue
    time.sleep(0.4)
print("done")
