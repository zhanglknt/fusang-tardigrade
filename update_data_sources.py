"""Rebuild DATA_SOURCES.json v2.0 with all new experiments."""
import json

data = {
    "_description": "Claim-to-data mapping for manuscript-auditor verification.",
    "_version": "2.0",
    "_updated": "2026-06-27: Added L3 validation, Mash benchmark, E2E classifier, Cohen's d correction",
    "_nrf_formula": "nRF = RF / (2*(n-3))",
    "_n200_max_rf": 394,
    "_cohens_d_formula": "d = (mean_x - mean_y) / pooled_sd",

    "tables": {
        "Table 1": {
            "description": "IMMI L0-1 vs MSA+ML across scales (30 seeds per condition)",
            "entries": [
                {"condition": "n=200 Clean", "immi_nrf": 0.102, "ft2_nrf": 0.096, "file": "benchmark_n200_clean_30seeds.csv"},
                {"condition": "n=200 Indel(0.02)", "immi_nrf": 0.078, "ft2_nrf": 0.080, "file": "benchmark_n200_indel_30seeds.csv"},
                {"condition": "n=500 Clean", "immi_nrf": 0.119, "ft2_nrf": 0.093, "file": "benchmark_n500_clean_30seeds.csv"},
                {"condition": "n=500 Indel", "immi_nrf": 0.095, "ft2_nrf": 0.083, "file": "benchmark_n500_indel_30seeds.csv"},
                {"condition": "n=1000 Clean", "immi_nrf": 0.115, "ft2_nrf": 0.091, "file": "benchmark_n1000_clean_30seeds.csv"}
            ]
        },
        "Table 2 (Indel scan)": {
            "description": "n=200, 130 seeds, MASTER_REBUILT source",
            "entries": [
                {"indel_rate": 0.005, "immi_nrf": 0.137},
                {"indel_rate": 0.010, "immi_nrf": 0.107},
                {"indel_rate": 0.020, "immi_nrf": 0.080},
                {"indel_rate": 0.050, "immi_nrf": 0.066}
            ],
            "file": "table3_corrected.csv"
        },
        "Table 3 (IQ-TREE2)": {
            "description": "IQ-TREE2 GTR vs k-mer methods",
            "entries": [
                {"method": "IMMI L0-1", "nrf": 0.080, "n_valid": 112},
                {"method": "FastTree2", "nrf": 0.085, "n_valid": 112},
                {"method": "IQ-TREE2 GTR", "nrf": 0.147, "n_valid": 121}
            ],
            "file": "iqtree2_benchmark_indel_n200.csv"
        },
        "Table 7 (Multi-k)": {
            "description": "Multi-k ensemble, n=200 indel, 30 seeds",
            "entries": [
                {"method": "k=5,gap2 spaced", "nrf": 0.112},
                {"method": "k=5 contiguous", "nrf": 0.099},
                {"method": "k=7 contiguous", "nrf": 0.102},
                {"method": "k=9 contiguous", "nrf": 0.112},
                {"method": "multi-k ensemble", "nrf": 0.105}
            ],
            "file": "benchmark_multik_ensemble_table7_full.csv"
        },
        "Table 8 (Competitors)": {
            "description": "n=200 indel=0.02, 27 seeds. Cohen's d CORRECTED: 8.6->20.15",
            "cohens_d_correction": {
                "claimed": 8.6,
                "corrected_TRUE_ref": 20.15,
                "corrected_FT2_ref": 1.99,
                "source_TRUE": "table8_results_TRUE_reference.csv",
                "note": "Original d=8.6 unreproducible from any data source"
            },
            "entries": [
                {"method": "IMMI L0-1", "nrf": 0.112},
                {"method": "Co-phylog", "nrf": 0.419},
                {"method": "KmerCosine k=5", "nrf": 0.099},
                {"method": "KmerCosine k=7", "nrf": 0.102}
            ]
        }
    },

    "supplements": {
        "L3_Validation": {
            "description": "L0 vs L1 vs L3/FT2 on n=200 indel, vs TRUE tree",
            "key_finding": "L1 (multi-k) = FT2 (MSA+ML), p<0.0001 vs L0",
            "entries": [
                {"method": "L0 (k=5 NJ)", "nrf_mean": 0.7428, "nrf_std": 0.0460, "n": 30},
                {"method": "L1 (multi-k NJ)", "nrf_mean": 0.5826, "nrf_std": 0.0445, "n": 30},
                {"method": "L3/FT2 (MAFFT+FastTree2)", "nrf_mean": 0.5919, "nrf_std": 0.0407, "n": 5, "note": "MAFFT flaky on Windows"}
            ],
            "file": "l3_validation_n200/l3_validation_results.json"
        },
        "Mash_Benchmark": {
            "description": "IMMI vs Mash on n=200, vs TRUE tree",
            "key_finding": "Mash nRF=1.005 on indels; spaced k-mer 1.35x better",
            "entries": [
                {"method": "IMMI Clean", "nrf_mean": 0.3760, "nrf_std": 0.0451, "n": 30},
                {"method": "IMMI Indel", "nrf_mean": 0.7418, "nrf_std": 0.0449, "n": 30},
                {"method": "Mash Clean", "nrf": 0.1624, "n": 1},
                {"method": "Mash Indel", "nrf": 1.0051, "n": 1}
            ],
            "file": "mash_vs_immi_results.json"
        },
        "E2E_Classifier": {
            "description": "Boundary classifier RF V4b, 88 scenarios",
            "key_finding": "88/88 correct, Wilson CI [0.958, 1.0]",
            "entries": [
                {"type": "coalescent", "expected": "STOP", "correct": "13/13"},
                {"type": "structured", "expected": "SPLIT", "correct": "75/75"},
                {"overall": "88/88, CI [0.958, 1.0]"}
            ],
            "file": "e2e_extended_results.json"
        }
    },

    "text_claims": [
        {"section": "Cohens d correction", "claim": "d=20.15 (not 8.6) for Co-phylog vs Fusang", "file": "table8_results_TRUE_reference.csv"},
        {"section": "Multi-k matches gold standard", "claim": "L1=0.583 vs FT2=0.592, multi-k NJ = MSA+ML", "file": "l3_validation_n200/l3_validation_results.json"},
        {"section": "Mash fails on indels", "claim": "Mash nRF=1.005 vs IMMI=0.742, spaced k-mer 1.35x better", "file": "mash_vs_immi_results.json"},
        {"section": "E2E classifier validation", "claim": "88/88 correct, CI [0.958, 1.0]", "file": "e2e_extended_results.json"},
        {"section": "IMMI L0-1 accuracy", "claim": "nRF=0.080 on n=200 indel", "file": "indel_benchmark_MASTER_REBUILT.csv"},
        {"section": "Indel robustness", "claim": "5.9% advantage at indel=0.02", "file": "table3_corrected.csv"},
        {"section": "IQ-TREE2 degradation", "claim": "1.8x worse than IMMI, d=3.1", "file": "iqtree2_benchmark_indel_n200.csv"}
    ],

    "_deprecated_files": [
        "benchmark_table8_ft2_ref.csv",
        "benchmark_table8_verification.csv",
        "benchmark_table8_n27_indel.csv",
        "benchmark_competitors_n200_indel.csv"
    ],
    "_deprecated_reason": "Old normalization max_rf=len(total), corrected to 2*(n-3)"
}

out_path = "repro_package/data/DATA_SOURCES.json"
with open(out_path, "w", encoding="utf-8") as f:
    json.dump(data, f, indent=2, ensure_ascii=False)

# Validate
with open(out_path, encoding="utf-8") as f:
    json.load(f)

print(f"DATA_SOURCES.json v2.0 written: {out_path}")
print(f"  Tables: {len(data['tables'])}")
print(f"  Supplements: {len(data['supplements'])}")
print(f"  Text claims: {len(data['text_claims'])}")
print("  JSON validation: OK")
