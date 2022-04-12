# WG3-pipeline-QTL
Part of the sc-eQTLgen consortium pipeline. Step 3, where the QTL mapping is done.
See the [wiki](https://github.com/sc-eQTLgen-consortium/WG3-pipeline-QTL/wiki) for more detailed information.

##General overview of the steps.

1) perform harmonization over all cells to cluster cell types (not sure how this overlaps with WG2)
2) determine which covariates to include in correction, per celltype, using harmonized dataset (e.g. PCs, sequencing-based covariates, sample quality measures, etc); not sure if this needs to be done over a combined dataset? we shortly discussed this last week...
3) perform an eQTL analysis within each dataset, using LIMIX pipeline (https://github.com/single-cell-genetics/limix_qtl)
4) perform 1000 permutations, summarize permutations per gene using the two beta-distribution parameters
5) collect summary statistics for step 3 and 4 at central site (per cell type and per dataset maximally 300 million cis-eQTLs sumtats + 40,000 permutation summary statistics, when using 20,000 genes; more likely half of this in single nucleus data)
6) meta-analyze non-permuted summary statistics per gene, combine beta distribution parameters from all datasets to adjust top non-permuted p-value per gene
7) determine false discovery rate using adjusted p-values per gene
