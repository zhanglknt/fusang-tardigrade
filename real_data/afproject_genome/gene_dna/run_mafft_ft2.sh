#!/bin/bash
# Step 3b: per-gene MAFFT alignment + FastTree2 ML tree
set -u
BASE="/mnt/d/系统发育树项目/Fusang/Fusang-main/real_data/afproject_genome/gene_dna"
mkdir -p "$BASE/aln" "$BASE/trees"
GENES="ND1 ND2 ND3 ND4 ND4L ND5 ND6 ATP6 ATP8 COX1 COX2 COX3 CYTB"
ENVBIN="/opt/miniconda3/envs/c_layer_ortho/bin"
for g in $GENES; do
  "$ENVBIN/mafft" --auto --quiet "$BASE/fasta/$g.fasta" > "$BASE/aln/$g.fasta" 2>/dev/null
  "$ENVBIN/FastTree" -nt -gtr -nosupport "$BASE/aln/$g.fasta" > "$BASE/trees/${g}_ft2.nwk" 2>/dev/null
  echo "$g: aln $(grep -c '>' "$BASE/aln/$g.fasta") seqs, tree $(wc -c < "$BASE/trees/${g}_ft2.nwk") bytes"
done
