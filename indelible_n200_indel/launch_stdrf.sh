#!/bin/bash
# Detached launcher: starts the standard-RF competitor rerun inside WSL
# fusang_bench with setsid+nohup so it survives the launching shell /
# agent session, then returns immediately.
#   progress : competitors_stdrf.log
#   resume   : .cache_competitors_n200_indel_130seeds_stdrf/  (per-seed JSON)
#   trees    : trees_stdrf/seed{N}_{method}.nwk
set -e
DIR="/mnt/d/系统发育树项目/Fusang/Fusang-main/indelible_n200_indel"
cd "$DIR"
setsid nohup bash run_competitors_stdrf.sh \
    > competitors_stdrf.log 2>&1 < /dev/null &
echo "launched pid=$!"
sleep 3
tail -5 competitors_stdrf.log 2>/dev/null || echo "(log not written yet)"
ps aux | grep "[r]un_competitors_stdrf" | head -5
