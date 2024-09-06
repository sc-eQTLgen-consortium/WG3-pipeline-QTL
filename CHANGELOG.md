# Changelog

All notable changes implemented in this branch compared to the main branch of [sc-eQTLgen-consortium](https://github.com/sc-eQTLgen-consortium/WG3-pipeline-QTL) (2023-08-22, SHA: `4a13d9386f70583dc2bde15477cd313659f48eba`) are documented in this file. 

Note that this branch is in beta and version 2.0.0 is not yet ready for release.


## [2.0.0] - eQTL mapping - 2023-09-21

#### Additions
- Dynamic time usage per rule
- Add settings to turn on / off certain calculations
- Added optional rule `filter_vcf` to filter the input VCF on samples you will be using for eQTL analysis (reduces memory and storage usage)
- Split `derive_expression_matrices` into `create_expression_matrices` and `prepare_eqtl_input` in order to be able to process all samples, even those missing genotype information
- Added example `wg2_pairing.csv` and `qc_mad_final.tab` files to this repo
- Moved `QC_functions.R` script to this repo
- Added merging of the pools and related rule `barcode_qc_plots` from WG1
- Added rules for sample eQTL QC (`kinship_pca`, `qtl_sample_qc`)
- Refactor rule `top_feature` into rule `minimal_postprocess` to allow output of all summary statistics
- Added printing of input settings for `run_qtl_mapping` rule for debugging
- Added `combine_qtl_mapping_runtime` rule for debugging

#### Fixes
- Fixed issue `processes_output.R` where `data.table::fread` fails to load complete gzipped files if too little tmp storage is available (see `Rdatatable` issue [#5095](https://github.com/Rdatatable/data.table/issues/5095)). This functionality has been reimplemented in python as `filter_variants.py`.

#### Changes
- Merged `Create_wp3_inputs` and `Qtl_Snakefile` into one snakemake
- Refactor code to (mostly) PEP8
- Improved clarity of rules by making input / output of rules more explicit
- Moved scripts outside of image for flexibility
- Moved all filepaths to settings file
- Moved all settings to settings file
- Made `ancestry` and `cell level` dynamic
- Refactor to enable parallel preprocessing of pools
- Refactor code to create seurat objects using `create_seurat` rule with all metadata files from WG1 and WG2 in a memory efficient way
- Made `kinship` optional as this functionality is now moved to WG1
- Refactor code to be more explicit on the LD / eQTL chunks it expects
- Split `derived_feature_annotation_and_chunks` into `feature_annotation` and `eqtl_chunks` in order to split eQTL from LD rules
- Refactor to redo PFlog1pPF as well as recalculate expression PCs after sample removal
- Made rule `compress_qtl` and `compress_ld` optional
