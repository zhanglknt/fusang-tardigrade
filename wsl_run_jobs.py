#!/usr/bin/env python3
"""
WSL-side batch driver for INDELible.

Runs inside Ubuntu-24.04 (INDELible is not installed on the Windows side).
Reads a job list produced by ``plan_indelible_jobs.py``, runs
``gen_indelible_sim.py`` for each job, and writes a per-job status file so
failures are recorded instead of being silently swallowed.

Usage (inside WSL):
    python3 /root/wsl_run_jobs.py /root/jobs.json /root/indel_bc_pool \\
            /root/job_status.json [--limit N] [--offset N]
"""

import argparse
import json
import os
import subprocess
import sys
import time

GEN = "/root/gen_indelible_sim.py"
WORK = "/root/indel_run_bc"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jobs")
    ap.add_argument("outdir")
    ap.add_argument("status")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--offset", type=int, default=0)
    args = ap.parse_args()

    with open(args.jobs) as fh:
        jobs = json.load(fh)
    jobs = jobs[args.offset:]
    if args.limit:
        jobs = jobs[:args.limit]

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)

    status = []
    t_start = time.time()
    for i, j in enumerate(jobs):
        cmd = [
            sys.executable, GEN,
            "--n", str(j["n"]),
            "--L", str(j["L"]),
            "--sub", repr(j["sub"]),
            "--indel", repr(j["indel"]),
            "--sub_scale", repr(j["sub_scale"]),
            "--seed", str(j["seed"]),
            "--outdir", args.outdir,
            "--workdir", WORK,
            "--tag", j["tag"],
            "--quiet",
        ]
        t0 = time.time()
        r = subprocess.run(cmd, capture_output=True, text=True)
        ok = (r.returncode == 0)
        rec = dict(j)
        rec["ok"] = ok
        rec["secs"] = round(time.time() - t0, 3)
        rec["err"] = "" if ok else ((r.stderr or "") + (r.stdout or ""))[-600:]
        status.append(rec)
        if not ok:
            sys.stderr.write("FAILED %s: %s\n" % (j["tag"], rec["err"][:300]))
        if (i + 1) % 100 == 0:
            el = time.time() - t_start
            sys.stderr.write("  %d/%d  %.0fs elapsed\n" % (i + 1, len(jobs), el))

    with open(args.status, "w") as fh:
        json.dump(status, fh)

    n_ok = sum(1 for s in status if s["ok"])
    print("DONE %d/%d ok in %.0fs" % (n_ok, len(status), time.time() - t_start))


if __name__ == "__main__":
    main()
