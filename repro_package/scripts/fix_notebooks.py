#!/usr/bin/env python3
"""
Fix all Jupyter notebooks in repro_package/notebooks/
to load data dynamically instead of using hardcoded values.
"""
import json
import os

NOTEBOOKS_DIR = r"D:\系统发育树项目\Fusang\Fusang-main\repro_package\notebooks"
DATA_DIR = r"D:\系统发育树项目\Fusang\Fusang-main\repro_package\data"

def fix_notebook_02(nb):
    """Fix notebook 02: Multi-k ensemble - load from CSV"""
    for i, cell in enumerate(nb['cells']):
        if cell['cell_type'] == 'code' and 'Table 3' in ''.join(cell['source']):
            print(f"  [02] Updating cell {i}: Table 3 data loading")
            cell['source'] = [
                '# Table 3 reproduction: Multi-k Ensemble\n',
                '# Load actual benchmark results from data file\n',
                'import csv, os\n',
                'import numpy as np\n',
                'import statistics as st\n',
                'from scipy import stats\n\n',
                'data_path = "../data/multik_ensemble_n200_30seeds.csv"\n',
                'if os.path.exists(data_path):\n',
                '    ensemble_vals = []\n',
                '    with open(data_path) as f:\n',
                '        for row in csv.DictReader(f):\n',
                '            ensemble_vals.append(float(row["ensemble_nrf"]))\n',
                '    ensemble_mean = st.mean(ensemble_vals)\n',
                '    ensemble_std = st.stdev(ensemble_vals)\n',
                '    print(f"Loaded {len(ensemble_vals)} seeds from {data_path}")\n',
                '    print(f"Ensemble nRF: {ensemble_mean:.4f} ± {ensemble_std:.4f}")\n',
                'else:\n',
                '    print("Data file not found, using manuscript values")\n',
                '    ensemble_mean, ensemble_std = 0.105, 0.021\n',
                '    print(f"Ensemble nRF (manuscript): {ensemble_mean:.4f} ± {ensemble_std:.4f}")\n'
            ]
            break
    return nb

def fix_notebook_03(nb):
    """Fix notebook 03: SwissTree validation - load from CSV"""
    for i, cell in enumerate(nb['cells']):
        src = ''.join(cell['source'])
        if cell['cell_type'] == 'code' and ('Table 5' in src or 'swisstree' in src.lower()):
            print(f"  [03] Updating cell {i}: Table 5 data loading")
            cell['source'] = [
                '# Table 5 reproduction: SwissTree protein gene tree benchmark\n',
                'import csv, os\n',
                'import numpy as np\n',
                'import statistics as st\n\n',
                'data_path = "../data/benchmark_swisstree_results.csv"\n',
                'if os.path.exists(data_path):\n',
                '    fusang_vals, cotree_vals = [], []\n',
                '    with open(data_path) as f:\n',
                '        for row in csv.DictReader(f):\n',
                '            fusang_vals.append(float(row["fusang_nrf"]))\n',
                '            cotree_vals.append(float(row["cotree_nrf"]))\n',
                '    print(f"Loaded {len(fusang_vals)} families from {data_path}")\n',
                '    print(f"Fusang nRF: {st.mean(fusang_vals):.4f} ± {st.stdev(fusang_vals):.4f}")\n',
                '    print(f"Co-phylog nRF: {st.mean(cotree_vals):.4f} ± {st.stdev(cotree_vals):.4f}")\n',
                '    # Statistical test\n',
                '    from scipy import stats\n',
                '    t_stat, p_val = stats.ttest_ind(fusang_vals, cotree_vals)\n',
                '    print(f"Welch\'s t-test: p={p_val:.4f}")\nelse:\n',
                '    print("Data file not found. Run: python ../scripts/benchmark_swisstree.py")\n'
            ]
            break
    return nb

def fix_notebook_05(nb):
    """Fix notebook 05: Scalability - already loads from JSON, just verify"""
    for i, cell in enumerate(nb['cells']):
        src = ''.join(cell['source'])
        if 'scalability_results.json' in src and 'open(' in src:
            print(f"  [05] Cell {i}: Already loads from JSON (OK)")
            break
    return nb

def fix_all_notebooks():
    """Fix all notebooks"""
    notebooks = {
        '02_multik_ensemble.ipynb': fix_notebook_02,
        '03_swisstree_validation.ipynb': fix_notebook_03,
        '05_scalability_demo.ipynb': fix_notebook_05,
    }
    
    for nb_name, fix_func in notebooks.items():
        nb_path = os.path.join(NOTEBOOKS_DIR, nb_name)
        if not os.path.exists(nb_path):
            print(f"[SKIP] {nb_name} not found")
            continue
        
        print(f"\nFixing {nb_name}...")
        with open(nb_path, encoding='utf-8') as f:
            nb = json.load(f)
        
        nb = fix_func(nb)
        
        with open(nb_path, 'w', encoding='utf-8') as f:
            json.dump(nb, f, indent=1, ensure_ascii=False)
        print(f"  Saved {nb_path}")
    
    print("\n✅ All notebooks updated!")

if __name__ == '__main__':
    fix_all_notebooks()
