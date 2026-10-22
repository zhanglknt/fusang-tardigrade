#!/bin/bash
# Recompute the competitor columns (INDELible n=200 indel, seeds 100-229) using
# the STANDARD root-independent Robinson-Foulds distance (nrf_standard.py)
# instead of the biased calc_nrf_simple.py convention.
#
# Pipelines / inputs / configurations are identical to the run that produced
# competitors_n200_indel_130seeds_fusangbench.csv; only the RF convention and
# the addition of an explicit cosine_k9 column change.  Every inferred tree is
# written to trees_stdrf/ so future re-scoring does not need a rerun.
#
#   bash run_competitors_stdrf.sh [seed_spec] [out_prefix] [workers]
set -e
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate fusang_bench
cd "/mnt/d/系统发育树项目/Fusang/Fusang-main"
SEEDS=${1:-100-229}
PREFIX=${2:-competitors_n200_indel_130seeds_stdrf}
WORKERS=${3:-2}
python -u run_competitors_stdrf.py --seeds "$SEEDS" --workers "$WORKERS" \
    --out-prefix "$PREFIX"
