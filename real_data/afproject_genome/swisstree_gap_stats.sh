#!/bin/bash
# Gap-content quantification: MAFFT --auto per SwissTree family, report gap fraction.
set -e
MAFFT=/opt/miniconda3/envs/c_layer_ortho/bin/mafft
SW='/mnt/d/系统发育树项目/Fusang/Fusang-main/real_data/swisstree/swisstree'
OUT='/mnt/d/系统发育树项目/Fusang/Fusang-main/real_data/afproject_genome/wsl_out'
mkdir -p ~/sw_gap && cd ~/sw_gap
rm -rf * 2>/dev/null || true
for fam in $(ls "$SW" | sed 's/_.*//' | sort -u); do
  rm -f ${fam}.fasta ${fam}.aln
  for f in "$SW"/${fam}_*.fasta; do
    awk -v id="$(basename "$f" .fasta)" '/^>/{print ">"id; next} {print}' "$f" >> ${fam}.fasta
  done
  $MAFFT --auto --quiet ${fam}.fasta > ${fam}.aln 2>/dev/null
  python3 - "$fam" << 'EOF'
import sys
fam = sys.argv[1]
seqs = []
for line in open(f'{fam}.aln'):
    if not line.startswith('>'):
        seqs.append(line.strip())
L = len(seqs[0]); gap = sum(s.count('-') for s in seqs)
print(f"{fam} n={len(seqs)} L={L} gap_fraction={gap/(L*len(seqs)):.4f}")
EOF
done
