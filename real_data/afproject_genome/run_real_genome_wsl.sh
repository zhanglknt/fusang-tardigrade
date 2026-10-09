#!/bin/bash
# External tools for real genome benchmarks (run inside WSL Ubuntu-24.04):
#   fish_mito: kmacs, mash, MAFFT+FastTree2 (MSA+ML), gap stats
#   ecoli_shigella: kmacs, mash (whole-genome MSA infeasible by design)
set -e
KMACS=~/af_tools/kmacs_src/kmacs
MASH=/opt/miniconda3/envs/skmer/bin/mash
MAFFT=/opt/miniconda3/envs/c_layer_ortho/bin/mafft
FASTTREE=/opt/miniconda3/envs/c_layer_ortho/bin/FastTree

SRC='/mnt/d/系统发育树项目/Fusang/Fusang-main/real_data/afproject_genome'
mkdir -p ~/real_bench && cd ~/real_bench
rm -rf * 2>/dev/null || true

# ---------- helper: merge per-genome fastas into one multi-fasta ----------
merge() {
  local dir=$1 out=$2
  rm -f "$out"
  for f in "$dir"/*.fasta; do
    awk -v id="$(basename "$f" .fasta)" '/^>/{print ">"id; next} {print}' "$f" >> "$out"
  done
}

# =========================== fish_mito ===========================
merge "$SRC/fish_mito" fish_all.fasta
echo "== fish_mito: kmacs =="
$KMACS -k 3 fish_all.fasta fish_kmacs_k3.dmat > /dev/null 2>&1
echo "== fish_mito: mash =="
$MASH sketch -o fish_msh -k 21 -s 10000 fish_all.fasta > /dev/null 2>&1
$MASH dist fish_msh.msh fish_msh.msh > fish_mash.tsv 2>/dev/null
echo "== fish_mito: MAFFT+FastTree2 =="
time $MAFFT --auto --quiet fish_all.fasta > fish_aln.fasta
$FASTTREE -nt -gtr -nosupport fish_aln.fasta > fish_ft2.nwk 2>/dev/null
# gap content of the alignment (excluding all-gap columns report too)
python3 - << 'EOF'
seqs, name = [], None
for line in open('fish_aln.fasta'):
    if line.startswith('>'): name = line[1:].strip()
    else: seqs.append(line.strip())
L = len(seqs[0]); gap = sum(s.count('-') for s in seqs)
allgap = sum(1 for c in zip(*seqs) if all(x == '-' for x in c))
print(f"FISH_ALIGNMENT n={len(seqs)} L={L} gap_fraction={gap/(L*len(seqs)):.4f} all_gap_cols={allgap}")
EOF

# ======================== ecoli_shigella ========================
merge "$SRC/ecoli_shigella" ecoli_all.fasta
echo "== ecoli_shigella: kmacs =="
time $KMACS -k 3 ecoli_all.fasta ecoli_kmacs_k3.dmat > /dev/null 2>&1
echo "== ecoli_shigella: mash =="
$MASH sketch -o ecoli_msh -k 21 -s 10000 ecoli_all.fasta > /dev/null 2>&1
$MASH dist ecoli_msh.msh ecoli_msh.msh > ecoli_mash.tsv 2>/dev/null

# ---- copy outputs back to Windows ----
mkdir -p "$SRC/wsl_out"
cp fish_kmacs_k3.dmat fish_mash.tsv fish_aln.fasta fish_ft2.nwk ecoli_kmacs_k3.dmat ecoli_mash.tsv "$SRC/wsl_out/"
echo "ALL WSL RUNS DONE"
