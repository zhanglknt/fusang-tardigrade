# Supplementary Table S17. Gene-length real DNA benchmark: 13 mitochondrial protein-coding genes from 25 fish genomes (AFproject fish_mito dataset)

Fusang: k-mer cosine distance + BioPython NJ (multi-k = mean of k=5,7,9 distance matrices). Anchor: MAFFT --auto + FastTree2 -nt -gtr -nosupport. Scoring: ETE3 Tree.compare (unrooted, max_rf = 2(n-3) = 44). Reference: Fischer et al. 2013 whole-mitogenome tree.

| Gene | Length mean bp (range) | Alignment gap frac | Multi-k vs ref | k=5 vs ref | k=7 vs ref | MAFFT+FT2 vs ref | Multi-k vs FT2 |
|------|----------------------|-------------------:|---------------:|-----------:|-----------:|-----------------:|---------------:|
| ND5 | 1840 (1837–1866) | 0.015 | 0.273 | 0.455 | 0.273 | 0.045 | 0.227 |
| COX1 | 1574 (1548–1596) | 0.014 | 0.409 | 0.591 | 0.409 | 0.182 | 0.364 |
| ND4 | 1381 (1380–1381) | 0.000 | 0.545 | 0.727 | 0.545 | 0.136 | 0.409 |
| CYTB | 1140 (1131–1142) | 0.001 | 0.591 | 0.636 | 0.636 | 0.364 | 0.500 |
| ND2 | 1046 (1045–1047) | 0.004 | 0.455 | 0.455 | 0.455 | 0.227 | 0.364 |
| ND1 | 975 (975–975) | 0.000 | 0.545 | 0.591 | 0.545 | 0.364 | 0.455 |
| COX3 | 785 (784–786) | 0.002 | 0.636 | 0.682 | 0.636 | 0.545 | 0.318 |
| COX2 | 691 (691–691) | 0.000 | 0.545 | 0.591 | 0.591 | 0.273 | 0.455 |
| ATP6 | 683 (683–684) | 0.001 | 0.591 | 0.636 | 0.545 | 0.318 | 0.545 |
| ND6 | 522 (519–531) | 0.017 | 0.500 | 0.591 | 0.318 | 0.273 | 0.455 |
| ND3 | 350 (349–352) | 0.012 | 0.682 | 0.727 | 0.727 | 0.591 | 0.727 |
| ND4L | 297 (297–297) | 0.000 | 0.500 | 0.409 | 0.545 | 0.364 | 0.455 |
| ATP8 | 168 (168–168) | 0.000 | 0.682 | 0.727 | 0.682 | 0.636 | 0.591 |
| **mean ± SD** | | | **0.535 ± 0.113** | 0.601 ± 0.107 | 0.531 ± 0.136 | **0.332 ± 0.175** | **0.451 ± 0.126** |

Notes: Fusang-vs-FT2 agreement (last column) is the primary methodological evidence, free of gene-tree/species-tree discordance; the discordance floor is set by the anchor itself (MAFFT+FT2 vs ref = 0.332 ± 0.175). Per-locus gap fractions 0–0.017 indicate a low-indel regime. Source: real_data/afproject_genome/gene_dna/gene_dna_results.json.
