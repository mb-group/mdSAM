# mdSAM: Discovery of self-assembly motifs regulating membrane-protein organization and function

Analysis notebooks and figure panels for the manuscript *Discovery of self-assembly motifs regulating membrane-protein organization and function* (Elazar, Trivedi, et al., Babu lab, St. Jude Children's Research Hospital).

## Abstract

While cell-surface clustering of membrane proteins is central to numerous processes, its molecular determinants remain elusive. We postulated that peptide motifs within disordered regions of transmembrane proteins may drive self-association and clustering and developed an approach to identify them. Analysis of the human proteome revealed that one-fifth of membrane proteins (945 proteins; >1,400 motifs), including receptors and transporters, contain motifs that we term mdSAMs (membrane-disordered Self-Assembly Motifs). Biophysical, biochemical, structural, and cell-based analyses of mdSAMs in model and native membrane proteins demonstrate that mdSAMs can drive self-assembly and clustering. In the serotonin receptor HTR7, a single mdSAM-disrupting substitution alters nanoscale organization, attenuates basal and agonist-induced signaling, and shifts the transcriptional state of the cell. Notably, mdSAM sequences are unique, which might enable discrimination of self from non‑self during self‑assembly, yet conserved among orthologs, underscoring their relevance. These findings reveal a widespread mechanism linking membrane-protein self-assembly to cell-surface organization and function.

## What is an mdSAM?

An **mdSAM** (membrane-disordered Self-Assembly Motif) is a short peptide motif, usually 7-9 residues (median 8), located in an intrinsically disordered, solvent-exposed region of a transmembrane protein, that is predicted to self-assemble through beta-strand-mediated, amyloid-like interactions. The working hypothesis is that mdSAMs provide weak, sequence-specific self-association that can cluster membrane proteins at the cell surface.

**How mdSAMs are identified**

1. A random-forest classifier is trained on WALTZ-DB 2.0 (1,402 hexapeptides: 510 amyloid-forming, 892 non-forming) using the one-hot-encoded hexapeptide sequence plus the summed beta-strand free energy (F1 = 0.77, stratified 5-fold cross-validation).
2. Every hexapeptide window in human membrane proteins is scored; positive windows are Self-Assembly Motifs (SAMs).
3. A SAM is called an mdSAM if it lies in a disordered region: AlphaFold pLDDT < 70, non-helical and solvent exposed, excluding transmembrane segments and signal peptides. A ROC-derived threshold on consecutive positive windows sets the high-confidence set (Fig. S1I).

**What was found.** 1,425 high-confidence mdSAMs in 945 membrane proteins (about 20% of the 4,886 human membrane proteins analysed), 91% of them in N- or C-terminal tails, and 72% of mdSAM-containing proteins carry a single motif. Of 44 synthesised mdSAM peptides, 36 (81%) gave Thioflavin-T signal above buffer; in the eight peptides tested, proline substitutions reduced the signal, abolishing it in seven; cryo-EM showed ordered fibrils. In cells, mdSAMs increase surface puncta density of a synthetic scaffold and of native full-length proteins. In the serotonin receptor HTR7, the V451P substitution in the C-terminal mdSAM (VSVWVLQS) changes nanoscale cluster organisation, lowers basal and agonist-induced cAMP signalling (Emax reduced by 25%) and shifts the transcriptional state of the cell. mdSAM sequences are diverse and unique, and conserved among orthologs.

> Naming: older notebooks call these motifs **dSAM** (disordered SAM) and refer to the structured-region counterpart as **mSAM/sSAM/msSAM**. dSAM in a file name or variable means mdSAM.

## Repository contents

```
notebooks/   the 23 analysis notebooks listed below, numbered in pipeline order (outputs cleared)
figures/     panel images exported from notebook outputs (one PNG per notebook cell)
FIGURE_MAP.csv   machine-readable map: paper panel -> notebook -> cell -> image
```

## What each notebook does

Cell numbers are 0-based indices into the notebook (as in the `cells` list of the `.ipynb`).

### Data preparation and classifier

| Notebook | Task | Paper panels | Original file name |
|---|---|---|---|
| `01_train_waltzdb_classifier.ipynb` | Trains and evaluates the random-forest classifier on WALTZ-DB 2.0 hexapeptides (one-hot sequence + summed beta-strand energy); reports test-set counts, stratified K-fold ROC and the beta-strand energy distributions. Saves the trained model. | Fig 1D; S1E-G | `class_waltzdb_for_paper.ipynb` |
| `02_alphafold_residue_table_and_early_plots.ipynb` | Builds the shared AlphaFold-derived residue table (`AF2_with_ENSG.pkl`, cell 205) used by most other notebooks; remaining cells hold earlier versions of exploratory plots. | (input table) | `Redo_Figures.ipynb` |
| `03_mdsam_sequence_logos.ipynb` | Sequence logos and composition of mdSAMs and SAMs in membrane proteins (logomaker); also builds `hTMs_uniprot_v1.pkl` and `AF_IDR_v3.pkl`. Contains additional gnomAD analyses that are not part of the paper figures. | S1H | `Seqlogo.ipynb` |

### Discovery and diversity (Figs 1-2, S1-S2)

| Notebook | Task | Paper panels | Original file name |
|---|---|---|---|
| `04_pasta_aggregation_comparison.ipynb` | Compares mdSAM predictions with PASTA energies and an aggregation threshold on windows around WALTZ-DB hexapeptides; ROC for the number of consecutive positive windows. | S1I-J | `PASTA_nHS_mdSAMs.ipynb` |
| `05_uniprot_keyword_enrichment.ipynb` | UniProt keyword enrichment (observed minus expected, FDR) for mdSAM-containing membrane proteins: biological process, molecular function, disease. | Fig 1G; S1K | `KeyWords_clustermap.ipynb` |
| `06_mdsam_sequence_space_embedding.ipynb` | Embeds the mdSAM sequence space from physicochemical features, colours it by composition, clusters it (GMM) with keyword enrichment, draws the charge-archetype ternary plot, the dot matrix of mdSAM-containing proteins and the pairwise Levenshtein-distance analysis. | Fig 1F, 2F-I, 3A; S2E-G | `umap_paper-Copy2.ipynb` |
| `07_hyperlopit_organelle_enrichment.ipynb` | Re-analyses the hyperLOPIT subcellular fractionation (U2OS): fraction profiles, t-SNE and Gaussian-mixture organelle assignment, and mdSAM enrichment per organelle. | Fig 2A; S2A-C | `HPA_fractions.ipynb` |
| `08_mdsam_topology_and_membrane_distance.ipynb` | Topology of mdSAM-containing proteins: number of TM segments (bitopic vs polytopic), N-/C-terminus vs loop, distance to the nearest membrane leaflet, number and length of mdSAMs. | Fig 2B-F | `mdSAM_TM_location.ipynb` |
| `09_mdsams_per_protein.ipynb` | Histogram of mdSAMs per protein (1, 2, 3 or more). | Fig 2E | `Number_of_mdSaM_per_protein.ipynb` |

### In vitro self-assembly (Fig 3, S3)

| Notebook | Task | Paper panels | Original file name |
|---|---|---|---|
| `10_tht_peptide_screen_96well.ipynb` | Thioflavin-T screen of mdSAM peptides in 96-well plates: per-peptide intensities, spectra, wild-type vs proline-mutant plates, overlay on the charge-archetype ternary plot. | Fig 3A-B; S3A, S3C | `ThT_96_well.ipynb` |
| `11_tht_gpcr_peptides_proline_mutants.ipynb` | ThT of GPCR mdSAM peptides (HTR7, CCR4, DRD1) and their proline mutants, including HTR7 V451P. | Fig 4C, 5B; S3C | `ThT_dSAM_mutations.ipynb` |
| `12_amytracker_peptides.ipynb` | Amytracker 520 fluorescence (100 uM) of wild-type and proline-mutant peptides. | Fig 4C; S3D | `Amytracker_ploting.ipynb` |

### Cell-based clustering (Figs 4-5, S5-S6)

| Notebook | Task | Paper panels | Original file name |
|---|---|---|---|
| `13_tirf_puncta_poisson_glm.ipynb` | Poisson generalised linear model on CellProfiler cell and puncta tables: puncta per cell area and rate ratio (WT/mutant) per protein and cell line (HEK, U2OS), SNAPtm motif swaps, full-length APP and HTR7. | Fig 4C-D, 5E; S5C | `Puncta_revision_after_madan_after_rejecten.ipynb` |
| `14_flow_cytometry_total_vs_surface.ipynb` | Flow cytometry (FCS files): total and surface SNAP signal for wild-type vs mutant constructs by DNA amount, with percentile gating against mock cells. | S6A | `total_vs_surface.ipynb` |
| `15_biogrid_interactions_by_mdsam.ipynb` | BioGRID analysis of membrane-protein interactions by mdSAM status: self-interactions and plasma-membrane interaction networks with hypergeometric and binomial tests. Contains further interaction analyses that are not part of the paper figures. | S6B-C | `BioGRID_organelle.ipynb` |
| `16_dstorm_htr7_cluster_analysis.ipynb` | dSTORM cluster analysis: ThunderSTORM localisation tables, Voronoi thresholding against randomised points, HDBSCAN clustering, cluster size and density stratified by local density, bootstrap of WT-vs-V451P skewness (tail-heaviness) for HTR7. Also holds earlier SNAPtm/APP analyses that are not part of the paper figures. | Fig 5F-G | `dSTROM_papar-Copy1.ipynb` |

### HTR7 signalling and transcription (Fig 6)

| Notebook | Task | Paper panels | Original file name |
|---|---|---|---|
| `17_htr7_camp_emax_ec50_and_expression_plots.ipynb` | HTR7 cAMP Emax and EC50 comparison between WT and V451P; per-protein total- and surface-expression plots with a mixed model adjusted for total expression. | Fig 6D; S6A | `mdSAM_plotting_notebook_v2.ipynb` |
| `18_htr7_rnaseq_pca_volcano.ipynb` | RNA-seq PCA and volcano plots from DESeq2 tables (HTR7 WT, V451P, empty vector). | Fig 6E | `ywang34_plotings.ipynb` |
| `19_htr7_rnaseq_gsea.ipynb` | Re-runs differential expression (pyDESeq2) and gene-set enrichment (gseapy) for the HTR7 RNA-seq data; heatmaps and volcano plots. | supports Fig 6E | `GSEA_analysis.ipynb` |

### Evolution (Fig 7, S7)

| Notebook | Task | Paper panels | Original file name |
|---|---|---|---|
| `20_conservation_alignment_based.ipynb` | Alignment-based conservation (Kullback-Leibler, Shannon, Lockless; pycanal) of mdSAMs vs disordered and structured residues; retention of mdSAMs in one-to-one orthologs vs alignment depth. | Fig 7B; S7A-B | `evolution_paper_version.ipynb` |
| `21_conservation_esm_embedding.ipynb` | Alignment-free conservation from protein language-model embeddings for mdSAMs, random IDR segments, ELM motifs and structured residues. | Fig 7B | `Conservarion_embedding.ipynb` |
| `22_ortholog_peptide_tht.ipynb` | ThT of NGFR and TNFRSF13C ortholog mdSAM peptides (wild type and proline mutant) across species. | Fig 7D; S7C | `ThT_orthologs.ipynb` |
| `23_overlap_with_larks_wirs_finches.ipynb` | Overlap of mdSAM-containing membrane proteins with LARKS, WIRS and FINCHES-predicted fuzzy interactions (UpSet plot). | S7D | `LARKS_WIRS_hyder.ipynb` |

## Where the paper figures come from

Each row points to the notebook and cell that produced the panel and to the exported image. The images are the inline notebook outputs, so fonts, panel layout and annotations differ from the final published figures, which were assembled in a vector editor.

### Figure 1: Discovery of mdSAMs

| Panel | Content | Notebook (cell) | Image |
|---|---|---|---|
| 1D | RF classifier test-set performance (amyloid vs non-amyloid) | `01_train_waltzdb_classifier.ipynb` (cell 26) | [1D_classifier_c26.png](figures/1D_classifier_c26.png) |
| 1F | Dot matrix: membrane proteins with predicted mdSAMs | `06_mdsam_sequence_space_embedding.ipynb` (cell 14) | [1F_embedding_c14.png](figures/1F_embedding_c14.png) |
| 1G | UniProt keyword categories (donut) | `05_uniprot_keyword_enrichment.ipynb` (cell 32) | [1G_keywords_c32.png](figures/1G_keywords_c32.png) |

### Figure 2: Diversity of mdSAMs

| Panel | Content | Notebook (cell) | Image |
|---|---|---|---|
| 2A | mdSAM enrichment by organelle (hyperLOPIT) | `07_hyperlopit_organelle_enrichment.ipynb` (cell 79) | [2A_hyperLOPIT_c79.png](figures/2A_hyperLOPIT_c79.png) |
| 2A | mdSAM enrichment by organelle (hyperLOPIT) | `07_hyperlopit_organelle_enrichment.ipynb` (cell 87) | [2A_hyperLOPIT_c87.png](figures/2A_hyperLOPIT_c87.png) |
| 2A | mdSAM enrichment by organelle (hyperLOPIT) | `07_hyperlopit_organelle_enrichment.ipynb` (cell 89) | [2A_hyperLOPIT_c89.png](figures/2A_hyperLOPIT_c89.png) |
| 2B | Bitopic vs polytopic enrichment | `08_mdsam_topology_and_membrane_distance.ipynb` (cell 70) | [2B_topology_c70.png](figures/2B_topology_c70.png) |
| 2B | Bitopic vs polytopic enrichment | `08_mdsam_topology_and_membrane_distance.ipynb` (cell 73) | [2B_topology_c73.png](figures/2B_topology_c73.png) |
| 2C | N-/C-terminus vs loop | `08_mdsam_topology_and_membrane_distance.ipynb` (cell 146) | [2C_topology_c146.png](figures/2C_topology_c146.png) |
| 2C | N-/C-terminus vs loop | `08_mdsam_topology_and_membrane_distance.ipynb` (cell 148) | [2C_topology_c148.png](figures/2C_topology_c148.png) |
| 2C | N-/C-terminus vs loop | `08_mdsam_topology_and_membrane_distance.ipynb` (cell 168) | [2C_topology_c168.png](figures/2C_topology_c168.png) |
| 2D | Distance to nearest membrane leaflet | `08_mdsam_topology_and_membrane_distance.ipynb` (cell 45) | [2D_topology_c45.png](figures/2D_topology_c45.png) |
| 2D | Distance to nearest membrane leaflet | `08_mdsam_topology_and_membrane_distance.ipynb` (cell 55) | [2D_topology_c55.png](figures/2D_topology_c55.png) |
| 2E | mdSAMs per protein | `09_mdsams_per_protein.ipynb` (cell 2) | [2E_nmdSAM_c2.png](figures/2E_nmdSAM_c2.png) |
| 2F | mdSAM length distribution | `08_mdsam_topology_and_membrane_distance.ipynb` (cell 76) | [2F_topology_c76.png](figures/2F_topology_c76.png) |
| 2G | Pairwise Levenshtein distance distribution / matrix | `06_mdsam_sequence_space_embedding.ipynb` (cell 207) | [2G_embedding_c207.png](figures/2G_embedding_c207.png) |
| 2G | Pairwise Levenshtein distance distribution / matrix | `06_mdsam_sequence_space_embedding.ipynb` (cell 211) | [2G_embedding_c211.png](figures/2G_embedding_c211.png) |
| 2H-I | Physicochemical embedding coloured by composition | `06_mdsam_sequence_space_embedding.ipynb` (cell 28) | [2HI_embedding_c28.png](figures/2HI_embedding_c28.png) |
| 2H-I | Physicochemical embedding coloured by composition | `06_mdsam_sequence_space_embedding.ipynb` (cell 37) | [2HI_embedding_c37.png](figures/2HI_embedding_c37.png) |

### Figure 3: mdSAM peptides self-assemble

| Panel | Content | Notebook (cell) | Image |
|---|---|---|---|
| 3A | Charge archetypes (ternary plot) | `10_tht_peptide_screen_96well.ipynb` (cell 25) | [3A_ThT96_c25.png](figures/3A_ThT96_c25.png) |
| 3A | Charge archetypes (ternary plot) | `10_tht_peptide_screen_96well.ipynb` (cell 26) | [3A_ThT96_c26.png](figures/3A_ThT96_c26.png) |
| 3B | ThT fluorescence per peptide with buffer/MAPT lines | `10_tht_peptide_screen_96well.ipynb` (cell 54) | [3B_ThT96_c54.png](figures/3B_ThT96_c54.png) |
| 3B | ThT fluorescence per peptide with buffer/MAPT lines | `10_tht_peptide_screen_96well.ipynb` (cell 73) | [3B_ThT96_c73.png](figures/3B_ThT96_c73.png) |

### Figure 4: mdSAMs promote clustering in cells

| Panel | Content | Notebook (cell) | Image |
|---|---|---|---|
| 4C | Amytracker WT vs proline peptides | `12_amytracker_peptides.ipynb` (cell 3) | [4C_Amytracker_c3.png](figures/4C_Amytracker_c3.png) |
| 4C | SNAPtm motif swaps: puncta density, HTR7 and KCNE4 mdSAM | `13_tirf_puncta_poisson_glm.ipynb` (cell 42) | [4C_puncta_c42.png](figures/4C_puncta_c42.png) |
| 4C | SNAPtm motif swaps: puncta density, HTR7 and KCNE4 mdSAM | `13_tirf_puncta_poisson_glm.ipynb` (cell 43) | [4C_puncta_c43.png](figures/4C_puncta_c43.png) |
| 4C / 5B | ThT of HTR7 WT vs V451P peptide | `11_tht_gpcr_peptides_proline_mutants.ipynb` (cell 16) | [4C_5B_ThT_GPCR_c16.png](figures/4C_5B_ThT_GPCR_c16.png) |
| 4C / 5B | ThT of HTR7 WT vs V451P peptide | `11_tht_gpcr_peptides_proline_mutants.ipynb` (cell 17) | [4C_5B_ThT_GPCR_c17.png](figures/4C_5B_ThT_GPCR_c17.png) |
| 4D / S5C | Full-length protein puncta density, Poisson GLM, HEK vs U2OS | `13_tirf_puncta_poisson_glm.ipynb` (cell 12) | [4D_S5C_puncta_c12.png](figures/4D_S5C_puncta_c12.png) |
| 4D / S5C | Full-length protein puncta density, Poisson GLM, HEK vs U2OS | `13_tirf_puncta_poisson_glm.ipynb` (cell 15) | [4D_S5C_puncta_c15.png](figures/4D_S5C_puncta_c15.png) |
| 4D / S5C | Full-length protein puncta density, Poisson GLM, HEK vs U2OS | `13_tirf_puncta_poisson_glm.ipynb` (cell 15) | [4D_S5C_puncta_c15_1.png](figures/4D_S5C_puncta_c15_1.png) |
| 4D / S5C | Full-length protein puncta density, Poisson GLM, HEK vs U2OS | `13_tirf_puncta_poisson_glm.ipynb` (cell 19) | [4D_S5C_puncta_c19.png](figures/4D_S5C_puncta_c19.png) |
| 4D / S5C | Full-length protein puncta density, Poisson GLM, HEK vs U2OS | `13_tirf_puncta_poisson_glm.ipynb` (cell 22) | [4D_S5C_puncta_c22.png](figures/4D_S5C_puncta_c22.png) |
| 4D / S5C | Full-length protein puncta density, Poisson GLM, HEK vs U2OS | `13_tirf_puncta_poisson_glm.ipynb` (cell 22) | [4D_S5C_puncta_c22_1.png](figures/4D_S5C_puncta_c22_1.png) |

### Figure 5: HTR7 mdSAM and nanoscale organisation

| Panel | Content | Notebook (cell) | Image |
|---|---|---|---|
| 5B | HTR7 and other GPCR peptide ThT spectra | `11_tht_gpcr_peptides_proline_mutants.ipynb` (cell 32) | [5B_ThT_GPCR_c32.png](figures/5B_ThT_GPCR_c32.png) |
| 5E | HTR7 WT vs V451P TIRF puncta | `13_tirf_puncta_poisson_glm.ipynb` (cell 25) | [5E_puncta_c25.png](figures/5E_puncta_c25.png) |
| 5F | HTR7 dSTORM localisation map and Voronoi/cluster overlay | `16_dstorm_htr7_cluster_analysis.ipynb` (cell 409) | [5F_dSTORM_c409.png](figures/5F_dSTORM_c409.png) |
| 5F | HTR7 dSTORM localisation map and Voronoi/cluster overlay | `16_dstorm_htr7_cluster_analysis.ipynb` (cell 409) | [5F_dSTORM_c409_1.png](figures/5F_dSTORM_c409_1.png) |
| 5F | HTR7 dSTORM localisation map and Voronoi/cluster overlay | `16_dstorm_htr7_cluster_analysis.ipynb` (cell 410) | [5F_dSTORM_c410.png](figures/5F_dSTORM_c410.png) |
| 5F | HTR7 dSTORM localisation map and Voronoi/cluster overlay | `16_dstorm_htr7_cluster_analysis.ipynb` (cell 410) | [5F_dSTORM_c410_1.png](figures/5F_dSTORM_c410_1.png) |
| 5G | Cluster stoichiometry by local density; bootstrap skewness | `16_dstorm_htr7_cluster_analysis.ipynb` (cell 440) | [5G_dSTORM_c440.png](figures/5G_dSTORM_c440.png) |
| 5G | Cluster stoichiometry by local density; bootstrap skewness | `16_dstorm_htr7_cluster_analysis.ipynb` (cell 451) | [5G_dSTORM_c451.png](figures/5G_dSTORM_c451.png) |
| 5G | Cluster stoichiometry by local density; bootstrap skewness | `16_dstorm_htr7_cluster_analysis.ipynb` (cell 452) | [5G_dSTORM_c452.png](figures/5G_dSTORM_c452.png) |
| 5G | Cluster stoichiometry by local density; bootstrap skewness | `16_dstorm_htr7_cluster_analysis.ipynb` (cell 453) | [5G_dSTORM_c453.png](figures/5G_dSTORM_c453.png) |

### Figure 6: HTR7 mdSAM, signalling and transcription

| Panel | Content | Notebook (cell) | Image |
|---|---|---|---|
| 6D | HTR7 cAMP Emax and EC50 summary statistics | `17_htr7_camp_emax_ec50_and_expression_plots.ipynb` (cell 18) | [6D_plotting_v2_c18.png](figures/6D_plotting_v2_c18.png) |
| 6E | RNA-seq PCA and V451P vs WT volcano | `18_htr7_rnaseq_pca_volcano.ipynb` (cell 3) | [6E_RNAseq_c3.png](figures/6E_RNAseq_c3.png) |
| 6E | RNA-seq PCA and V451P vs WT volcano | `18_htr7_rnaseq_pca_volcano.ipynb` (cell 8) | [6E_RNAseq_c8.png](figures/6E_RNAseq_c8.png) |

### Figure 7: Evolutionary conservation

| Panel | Content | Notebook (cell) | Image |
|---|---|---|---|
| 7B | ESM2 embedding conservation by segment class | `21_conservation_esm_embedding.ipynb` (cell 15) | [7B_ESMcons_c15.png](figures/7B_ESMcons_c15.png) |
| 7B | Conservation by segment class (related analyses) | `20_conservation_alignment_based.ipynb` (cell 22) | [7B_MSAcons_c22.png](figures/7B_MSAcons_c22.png) |
| 7B | Conservation by segment class (related analyses) | `20_conservation_alignment_based.ipynb` (cell 28) | [7B_MSAcons_c28.png](figures/7B_MSAcons_c28.png) |
| 7D | Ortholog ThT: WT vs proline mutant per species | `22_ortholog_peptide_tht.ipynb` (cell 6) | [7D_ThT_orth_c6.png](figures/7D_ThT_orth_c6.png) |
| 7D | Ortholog ThT: WT vs proline mutant per species | `22_ortholog_peptide_tht.ipynb` (cell 16) | [7D_ThT_orth_c16.png](figures/7D_ThT_orth_c16.png) |

### Supplementary Figure S1: Validation of mdSAM prediction

| Panel | Content | Notebook (cell) | Image |
|---|---|---|---|
| S1F | ROC, stratified 5-fold CV | `01_train_waltzdb_classifier.ipynb` (cell 15) | [S1F_classifier_c15.png](figures/S1F_classifier_c15.png) |
| S1G | Beta-strand energetics | `01_train_waltzdb_classifier.ipynb` (cell 12) | [S1G_classifier_c12.png](figures/S1G_classifier_c12.png) |
| S1H | Sequence logos, mdSAMs vs SAMs | `03_mdsam_sequence_logos.ipynb` (cell 35) | [S1H_seqlogo_c35.png](figures/S1H_seqlogo_c35.png) |
| S1H | Sequence logos, mdSAMs vs SAMs | `03_mdsam_sequence_logos.ipynb` (cell 37) | [S1H_seqlogo_c37.png](figures/S1H_seqlogo_c37.png) |
| S1I | ROC for number of consecutive windows classified as amyloid | `04_pasta_aggregation_comparison.ipynb` (cell 17) | [S1I_PASTA_c17.png](figures/S1I_PASTA_c17.png) |
| S1I | ROC for number of consecutive windows classified as amyloid | `04_pasta_aggregation_comparison.ipynb` (cell 21) | [S1I_PASTA_c21.png](figures/S1I_PASTA_c21.png) |
| S1J | PASTA energy of mdSAMs | `04_pasta_aggregation_comparison.ipynb` (cell 14) | [S1J_PASTA_c14.png](figures/S1J_PASTA_c14.png) |
| S1K | UniProt keyword enrichment (process, function, disease) | `05_uniprot_keyword_enrichment.ipynb` (cell 47) | [S1K_keywords_c47.png](figures/S1K_keywords_c47.png) |
| S1K | UniProt keyword enrichment (process, function, disease) | `05_uniprot_keyword_enrichment.ipynb` (cell 47) | [S1K_keywords_c47_1.png](figures/S1K_keywords_c47_1.png) |
| S1K | UniProt keyword enrichment (process, function, disease) | `05_uniprot_keyword_enrichment.ipynb` (cell 73) | [S1K_keywords_c73.png](figures/S1K_keywords_c73.png) |
| S1K | UniProt keyword enrichment (process, function, disease) | `05_uniprot_keyword_enrichment.ipynb` (cell 74) | [S1K_keywords_c74.png](figures/S1K_keywords_c74.png) |

### Supplementary Figure S2: mdSAM diversity and localisation

| Panel | Content | Notebook (cell) | Image |
|---|---|---|---|
| S2A | hyperLOPIT profiles per compartment | `07_hyperlopit_organelle_enrichment.ipynb` (cell 26) | [S2A_hyperLOPIT_c26.png](figures/S2A_hyperLOPIT_c26.png) |
| S2A | hyperLOPIT profiles per compartment | `07_hyperlopit_organelle_enrichment.ipynb` (cell 29) | [S2A_hyperLOPIT_c29.png](figures/S2A_hyperLOPIT_c29.png) |
| S2B-C | hyperLOPIT t-SNE and GMM organelle assignment | `07_hyperlopit_organelle_enrichment.ipynb` (cell 48) | [S2BC_hyperLOPIT_c48.png](figures/S2BC_hyperLOPIT_c48.png) |
| S2B-C | hyperLOPIT t-SNE and GMM organelle assignment | `07_hyperlopit_organelle_enrichment.ipynb` (cell 68) | [S2BC_hyperLOPIT_c68.png](figures/S2BC_hyperLOPIT_c68.png) |
| S2E-F | Amino-acid composition and organelle on physicochemical map | `06_mdsam_sequence_space_embedding.ipynb` (cell 101) | [S2EF_embedding_c101.png](figures/S2EF_embedding_c101.png) |
| S2E-F | Amino-acid composition and organelle on physicochemical map | `06_mdsam_sequence_space_embedding.ipynb` (cell 122) | [S2EF_embedding_c122.png](figures/S2EF_embedding_c122.png) |
| S2E-F | Amino-acid composition and organelle on physicochemical map | `06_mdsam_sequence_space_embedding.ipynb` (cell 124) | [S2EF_embedding_c124.png](figures/S2EF_embedding_c124.png) |
| S2G | GMM groups and keyword enrichment | `06_mdsam_sequence_space_embedding.ipynb` (cell 49) | [S2G_embedding_c49.png](figures/S2G_embedding_c49.png) |
| S2G | GMM groups and keyword enrichment | `06_mdsam_sequence_space_embedding.ipynb` (cell 55) | [S2G_embedding_c55.png](figures/S2G_embedding_c55.png) |
| S2G | GMM groups and keyword enrichment | `06_mdsam_sequence_space_embedding.ipynb` (cell 65) | [S2G_embedding_c65.png](figures/S2G_embedding_c65.png) |
| S2G | GMM groups and keyword enrichment | `06_mdsam_sequence_space_embedding.ipynb` (cell 67) | [S2G_embedding_c67.png](figures/S2G_embedding_c67.png) |
| S2G | GMM groups and keyword enrichment | `06_mdsam_sequence_space_embedding.ipynb` (cell 69) | [S2G_embedding_c69.png](figures/S2G_embedding_c69.png) |

### Supplementary Figure S3: Peptide assembly

| Panel | Content | Notebook (cell) | Image |
|---|---|---|---|
| S3A | ThT spectra, all peptides | `10_tht_peptide_screen_96well.ipynb` (cell 86) | [S3A_ThT96_c86.png](figures/S3A_ThT96_c86.png) |
| S3C | ThT spectra, WT vs proline mutants | `10_tht_peptide_screen_96well.ipynb` (cell 80) | [S3C_ThT96_c80.png](figures/S3C_ThT96_c80.png) |
| S3D | Amytracker WT vs proline | `12_amytracker_peptides.ipynb` (cell 3) | [S3D_Amytracker_c3.png](figures/S3D_Amytracker_c3.png) |

### Supplementary Figure S6: Expression controls and BioGRID analysis

| Panel | Content | Notebook (cell) | Image |
|---|---|---|---|
| S6A | Total expression bar plots, per protein | `17_htr7_camp_emax_ec50_and_expression_plots.ipynb` (cell 13) | [S6A_plotting_v2_c13.png](figures/S6A_plotting_v2_c13.png) |
| S6A | Total expression bar plots, per protein | `17_htr7_camp_emax_ec50_and_expression_plots.ipynb` (cell 17) | [S6A_plotting_v2_c17.png](figures/S6A_plotting_v2_c17.png) |
| S6A | Total expression, WT vs proline, by DNA amount | `14_flow_cytometry_total_vs_surface.ipynb` (cell 45) | [S6A_total_vs_surface_c45.png](figures/S6A_total_vs_surface_c45.png) |
| S6A | Total expression, WT vs proline, by DNA amount | `14_flow_cytometry_total_vs_surface.ipynb` (cell 98) | [S6A_total_vs_surface_c98.png](figures/S6A_total_vs_surface_c98.png) |
| S6B | BioGRID self-interactions, mdSAM+ (n=658) and mdSAM- (n=2,642) | `15_biogrid_interactions_by_mdsam.ipynb` (cell 171) | [S6B_BioGRID_c171.png](figures/S6B_BioGRID_c171.png) |
| S6B | BioGRID self-interactions, mdSAM+ (n=658) and mdSAM- (n=2,642) | `15_biogrid_interactions_by_mdsam.ipynb` (cell 172) | [S6B_BioGRID_c172.png](figures/S6B_BioGRID_c172.png) |
| S6C | Plasma-membrane PPIs by mdSAM status | `15_biogrid_interactions_by_mdsam.ipynb` (cell 78) | [S6C_BioGRID_c78.png](figures/S6C_BioGRID_c78.png) |
| S6C | Plasma-membrane PPIs by mdSAM status | `15_biogrid_interactions_by_mdsam.ipynb` (cell 81) | [S6C_BioGRID_c81.png](figures/S6C_BioGRID_c81.png) |
| S6C | Plasma-membrane PPIs by mdSAM status | `15_biogrid_interactions_by_mdsam.ipynb` (cell 84) | [S6C_BioGRID_c84.png](figures/S6C_BioGRID_c84.png) |
| S6C | Plasma-membrane PPIs by mdSAM status | `15_biogrid_interactions_by_mdsam.ipynb` (cell 86) | [S6C_BioGRID_c86.png](figures/S6C_BioGRID_c86.png) |

### Supplementary Figure S7: Conservation

| Panel | Content | Notebook (cell) | Image |
|---|---|---|---|
| S7A | Alignment-based conservation metrics by residue class | `20_conservation_alignment_based.ipynb` (cell 24) | [S7A_MSAcons_c24.png](figures/S7A_MSAcons_c24.png) |
| S7A | Alignment-based conservation metrics by residue class | `20_conservation_alignment_based.ipynb` (cell 25) | [S7A_MSAcons_c25.png](figures/S7A_MSAcons_c25.png) |
| S7B | mdSAM retention in orthologs vs alignment depth | `20_conservation_alignment_based.ipynb` (cell 53) | [S7B_MSAcons_c53.png](figures/S7B_MSAcons_c53.png) |
| S7B | mdSAM retention in orthologs vs alignment depth | `20_conservation_alignment_based.ipynb` (cell 71) | [S7B_MSAcons_c71.png](figures/S7B_MSAcons_c71.png) |
| S7C | Ortholog ThT spectra | `22_ortholog_peptide_tht.ipynb` (cell 20) | [S7C_ThT_orth_c20.png](figures/S7C_ThT_orth_c20.png) |
| S7C | Ortholog ThT spectra | `22_ortholog_peptide_tht.ipynb` (cell 20) | [S7C_ThT_orth_c20_1.png](figures/S7C_ThT_orth_c20_1.png) |
| S7C | Ortholog ThT spectra | `22_ortholog_peptide_tht.ipynb` (cell 21) | [S7C_ThT_orth_c21.png](figures/S7C_ThT_orth_c21.png) |
| S7C | Ortholog ThT spectra | `22_ortholog_peptide_tht.ipynb` (cell 21) | [S7C_ThT_orth_c21_1.png](figures/S7C_ThT_orth_c21_1.png) |
| S7D | UpSet: overlap with LARKs, WIRS, FINCHES | `23_overlap_with_larks_wirs_finches.ipynb` (cell 16) | [S7D_UpSet_c16.png](figures/S7D_UpSet_c16.png) |

Panels that are schematics, microscopy or EM images, NMR/mass-photometry data, or were produced outside these notebooks are not listed.

## Running the notebooks

The notebooks are distributed with their outputs cleared; the panel images they produced are in `figures/`. The last column of the table above gives each notebook's original file name.

The notebooks were written for an interactive session and read intermediate files (for example `AF2_with_ENSG.pkl`, `hTMs_uniprot_v1.pkl`, `dSAM_MPs_v_paper_v2_new_locs.pkl`) that are not stored in this repository. Public inputs are WALTZ-DB 2.0, the AlphaFold human-proteome table of N. Davey (github.com/normandavey/ProcessedAlphafold), UniProt, the hyperLOPIT data (Thul et al., 2017), BioGRID 4.4.x, PANTHER 17 and GTEx v8. Supplementary data files accompany the paper.

Python packages imported by these notebooks (pip names can differ from import names): `Bio`, `Levenshtein`, `PIL`, `PyPDF2`, `adjustText`, `bokeh`, `cmasher`, `fcsparser`, `geopandas`, `gseapy`, `holoviews`, `intersect`, `iupred2a_lib`, `joblib`, `localcider`, `logomaker`, `matplotlib`, `matplotlib_venn`, `mycolorpy`, `networkx`, `numpy`, `pandas`, `plotly`, `pyaaisc`, `pycanal`, `pydeseq2`, `radialtree`, `rdfpy`, `ripleyk`, `sanbomics`, `scanpy`, `scipy`, `seaborn`, `shapely`, `sklearn`, `statsmodels`, `ternary`, `tqdm`, `umap`, `uncertainties`, `upsetplot`. `iupred2a_lib`, `pyaaisc` and `intersect` are local modules rather than PyPI packages.

## Citation and contact

Citation: to be added after publication. Correspondence: M. Madan Babu (madan.babu@stjude.org) and Assaf Elazar (assaf.elazar@stjude.org).
