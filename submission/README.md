NAR Submission Package
=====================
Fusang: Tardigrade Edition
Spaced k-mer Alignment-Free Phylogenetic Inference Resilient to Indel-Rich Sequence Evolution

Submission Date: June 27, 2026
Journal: Nucleic Acids Research (NAR) — Methods Article

---

## Package Contents

```
submission/
├── README.md                                  (this file)
├── NAR_Manuscript_Fusang_Tardigrade.docx       Main manuscript (MS Word)
├── Cover_Letter_NAR.docx                       Cover letter to editors
├── Supplementary_Material.docx                 Supplementary tables, figures, and notes
├── figures/
│   ├── main/                                   Figure 1-6 files (PNG)
│   └── supplementary/                          Figure S1-S9 files (PNG)
├── generate_cover_letter.js                    Script to regenerate cover letter
├── generate_supplementary.js                   Script to regenerate supplementary material
└── repro_package_v1.3.zip                      Full reproducibility package (4.8 MB)
```

## Key Findings Highlights

1. **Multi-k NJ ≈ MSA+ML on indels** (nRF=0.583 vs 0.592, 30-seed, vs TRUE tree)
2. **Mash collapses on indels** (nRF=1.005 = random) while spaced cosine degrades gracefully
3. **E2E classifier 100% accuracy** (88/88 scenarios, Wilson CI [0.958, 1.0])
4. **Cohen's d=20.15** (Co-phylog vs Fusang, TRUE reference)

## Reproducibility

All benchmark data, analysis scripts, and figure generation code are included in
repro_package_v1.3.zip. See DATA_SOURCES.json for complete data provenance mapping.

## Contact

Corresponding author: [to be filled]
GitHub: https://github.com/fusang-dev/fusang-tardigrade
