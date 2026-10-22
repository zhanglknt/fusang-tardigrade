#!/bin/bash
# Wait (max ~9.3 min) for the stdrf competitor run to finish, then print status.
cd "/mnt/d/系统发育树项目/Fusang/Fusang-main/indelible_n200_indel" || exit 9
for i in $(seq 1 55); do
    if ! pgrep -f "run_competitors_stdrf.py" > /dev/null; then
        echo "=== RUN FINISHED ==="
        break
    fi
    sleep 10
done
echo "seeds_done=$(grep -cE '^\[' competitors_stdrf.log)"
grep -E "^\[" competitors_stdrf.log | tail -2
if pgrep -f "run_competitors_stdrf.py" > /dev/null; then
    echo "still_running=yes"
else
    echo "still_running=no"
fi
