#!/usr/bin/env python3
"""Step 2: extract 13 mtDNA protein-coding genes per genome from GenBank."""
import os, sys, json
from Bio import SeqIO

BASE = os.path.dirname(os.path.abspath(__file__))
GBK = os.path.join(BASE, "gbk")
FASTA = os.path.join(BASE, "fasta")
os.makedirs(FASTA, exist_ok=True)

GENES = ["ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
         "ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB"]
SYN = {}
for g in GENES:
    SYN[g.lower()] = g
# common synonyms seen in mt genbank files
SYN.update({
    "nad1": "ND1", "nd1": "ND1", "nad2": "ND2", "nad3": "ND3",
    "nad4": "ND4", "nad4l": "ND4L", "nad5": "ND5", "nad6": "ND6",
    "atp6": "ATP6", "atpase6": "ATP6", "atp8": "ATP8", "atpase8": "ATP8",
    "cox1": "COX1", "co1": "COX1", "coi": "COX1", "coxi": "COX1",
    "cox2": "COX2", "co2": "COX2", "coii": "COX2", "coxii": "COX2",
    "cox3": "COX3", "co3": "COX3", "coiii": "COX3", "coxiii": "COX3",
    "cytb": "CYTB", "cob": "CYTB", "cytochrome b": "CYTB",
})

accs = [fn[:-4] for fn in sorted(os.listdir(GBK)) if fn.endswith(".gbk")]
per_gene = {g: {} for g in GENES}   # gene -> acc -> seq
missing = {g: [] for g in GENES}

for acc in accs:
    rec = SeqIO.read(os.path.join(GBK, acc + ".gbk"), "genbank")
    found = {}
    for feat in rec.features:
        if feat.type != "CDS":
            continue
        gname = None
        if "gene" in feat.qualifiers:
            gname = feat.qualifiers["gene"][0].strip().lower()
        if gname not in SYN and "product" in feat.qualifiers:
            # try to match product like "NADH dehydrogenase subunit 1"
            prod = feat.qualifiers["product"][0].strip().lower()
            prodmap = {
                "nadh dehydrogenase subunit 1": "ND1",
                "nadh dehydrogenase subunit 2": "ND2",
                "nadh dehydrogenase subunit 3": "ND3",
                "nadh dehydrogenase subunit 4": "ND4",
                "nadh dehydrogenase subunit 4l": "ND4L",
                "nadh dehydrogenase subunit 5": "ND5",
                "nadh dehydrogenase subunit 6": "ND6",
                "atp synthase f0 subunit 6": "ATP6", "atp synthase subunit 6": "ATP6",
                "atp synthase f0 subunit 8": "ATP8", "atp synthase subunit 8": "ATP8",
                "cytochrome c oxidase subunit i": "COX1",
                "cytochrome c oxidase subunit ii": "COX2",
                "cytochrome c oxidase subunit iii": "COX3",
                "cytochrome b": "CYTB",
            }
            gname = prodmap.get(prod)
        canon = SYN.get(gname) if gname else None
        if canon and canon not in found:
            seq = str(feat.extract(rec.seq))
            found[canon] = seq
    for g in GENES:
        if g in found:
            per_gene[g][acc] = found[g]
        else:
            missing[g].append(acc)

# write per-gene multi-fasta
stats = {}
for g in GENES:
    path = os.path.join(FASTA, g + ".fasta")
    lens = []
    with open(path, "w") as f:
        for acc in accs:
            if acc in per_gene[g]:
                s = per_gene[g][acc]
                lens.append(len(s))
                f.write(f">{acc}\n{s}\n")
    stats[g] = {"n": len(per_gene[g]), "len_min": min(lens) if lens else 0,
                "len_max": max(lens) if lens else 0,
                "len_mean": round(sum(lens)/len(lens), 1) if lens else 0,
                "missing": missing[g]}
    print(f"{g:6s} n={stats[g]['n']:2d} len {stats[g]['len_min']}-{stats[g]['len_max']}"
          + (f"  MISSING: {missing[g]}" if missing[g] else ""))

json.dump(stats, open(os.path.join(BASE, "extract_stats.json"), "w"), indent=1)
print("saved extract_stats.json")
