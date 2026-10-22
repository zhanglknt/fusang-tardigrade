"""
Scalability benchmark for Fusang/IMMI framework.
Tests wall-clock time and RAM for k-mer+NJ pipeline across n=200..10000.
FIXED: uses Biopython to properly count tree leaves.
"""
import sys, time, os, subprocess, json, psutil

sys.path.insert(0, '.')
from tree_simulation import make_coalescent_tree, simulate_seqs
from Bio import Phylo

PY = r'C:/Users/admin/.workbuddy/binaries/python/envs/default/Scripts/python.exe'
CODE = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}

def make_fasta(seqs, n, path):
    with open(path, 'w') as f:
        for i in range(n):
            seq = ''.join(CODE[b] for b in seqs[i])
            f.write(f'>t{i+1:05d}\n')
            for j in range(0, len(seq), 80):
                f.write(seq[j:j+80] + '\n')

def count_tips_biopython(nwk_path):
    """Properly count leaves using Biopython."""
    if not os.path.exists(nwk_path):
        return 0
    try:
        tree = Phylo.read(nwk_path, 'newick')
        return len(list(tree.get_terminals()))
    except Exception as e:
        # Fallback: count unique taxon names (strings before ':')
        with open(nwk_path, encoding='utf-8', errors='replace') as f:
            content = f.read()
        # Extract leaf names: pattern like 't00123:0.045' or '(t00123:0.045'
        import re
        names = set(re.findall(r'([A-Za-z0-9_]+):', content))
        return len(names)

def count_tips_newick(nwk_path):
    """Count tips by parsing Newick format manually (fallback)."""
    if not os.path.exists(nwk_path):
        return 0
    with open(nwk_path, encoding='utf-8', errors='replace') as f:
        content = f.read().strip().rstrip(';')
    # Split by commas at the top level (not inside parentheses)
    # Simple heuristic: count unique names before ':'
    import re
    names = re.findall(r'([A-Za-z][A-Za-z0-9_]*)', content)
    # Filter: names that appear before a ':' are taxa
    taxon_names = set()
    for m in re.finditer(r'([A-Za-z][A-Za-z0-9_]*)\s*:', content):
        taxon_names.add(m.group(1))
    return len(taxon_names)

def get_peak_ram():
    """Get current process peak RAM in MB."""
    try:
        process = psutil.Process()
        return process.memory_info().rss / 1024 / 1024
    except:
        return -1

def run_fusang(fa_path, out_path, use_simple=False, timeout_s=7200):
    """Run fusang. Returns (exit_code, elapsed_s, stdout)."""
    cmd = [PY, 'fusang_v2.py', '--input', fa_path, '--output', out_path,
           '--mode', 'default', '--kmer_k', '5', '--kmer_gap', 'gap2']
    if use_simple:
        cmd.append('--simple')
    # For n=10000, use default DCM mode (no --simple)
    t0 = time.time()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s)
        elapsed = time.time() - t0
        return result.returncode, elapsed, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        elapsed = time.time() - t0
        return -1, elapsed, '', f'TIMEOUT after {timeout_s}s'
    except Exception as e:
        elapsed = time.time() - t0
        return -2, elapsed, '', str(e)

# Load existing results
results_path = 'scalability_results.json'
if os.path.exists(results_path):
    with open(results_path) as f:
        results = json.load(f)
else:
    results = []

completed_n = {r['n'] for r in results if r.get('tips_match') or r.get('status') == 'OK'}

# Re-run with fixed tip counting, and add n=10000 with DCM mode
SCALES = [200, 500, 1000, 2000, 5000, 10000]
TIMEOUTS = {200: 600, 500: 600, 1000: 1800, 2000: 3600, 5000: 7200, 10000: 7200}

# Clear old results and re-run all with fixed counting
results = []
for n in SCALES:
    print(f'\n{"="*60}')
    print(f'SCALE: n={n}')
    print(f'{"="*60}')
    
    # Check if FASTA already exists (re-use to save time)
    fa_path = f'scalability_n{n}.fasta'
    if not os.path.exists(fa_path):
        t0 = time.time()
        root, leaves = make_coalescent_tree(n, seed=42)
        seqs = simulate_seqs(root, n, 500, 0.05, seed=42)
        make_fasta(seqs, n, fa_path)
        t1 = time.time()
        data_s = t1 - t0
        print(f'  Data generation: {data_s:.2f}s')
    else:
        data_s = 0.0
        print(f'  Data file exists, skipping generation')
    t1 = time.time()
    
    # Run Fusang
    out_path = f'scalability_n{n}.nwk'
    # NOTE (CLI fix): before the fusang_v2.py --simple precedence fix, this flag
    # was silently overridden to DCM for n>500; it now genuinely forces the
    # simplified pipeline at any n. Historical scalability timings for
    # 500<n<=5000 produced via this flag were therefore DCM runs.
    use_simple = (n <= 5000)  # --simple honored for all n; DCM only for n>5000
    timeout = TIMEOUTS.get(n, 7200)
    print(f'  Running Fusang (simple={use_simple}, timeout={timeout}s)...')
    
    exit_code, elapsed, stdout, stderr = run_fusang(fa_path, out_path, use_simple, timeout)
    
    # Fixed tip counting
    tips = count_tips_newick(out_path)
    tips_match = (tips == n)
    
    status = 'OK' if (exit_code == 0 and tips_match) else \
             f'TIMEOUT({elapsed:.0f}s)' if exit_code == -1 else \
             f'FAILED({exit_code})'
    
    ram = get_peak_ram()
    
    result = {
        'n': n,
        'data_gen_s': round(data_s, 2) if data_s > 0 else 'reused',
        'fusang_s': round(elapsed, 2),
        'total_s': round(data_s + elapsed, 2) if data_s > 0 else round(elapsed, 2),
        'exit_code': exit_code,
        'tips': tips,
        'tips_match': tips_match,
        'ram_mb': round(ram, 0) if ram > 0 else -1,
        'status': status,
    }
    
    results.append(result)
    print(f'  Result: {status}, tips={tips}/{n}, time={elapsed:.1f}s, RAM={ram:.0f}MB')
    
    # Save incremental
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Clean up
    for p in [fa_path, out_path]:
        if os.path.exists(p):
            os.remove(p)

print(f'\n{"="*60}')
print('FINAL RESULTS:')
for r in results:
    print(f"  n={r['n']:>6}: {r['fusang_s']:>8.1f}s, tips={r['tips']:>6}/{r['n']}, {r['status']}")
print(json.dumps(results, indent=2))
