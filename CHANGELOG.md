# Changelog

All notable changes implemented in this branch compared to the main branch of [sc-eQTLgen-consortium](https://github.com/sc-eQTLgen-consortium/WG3-pipeline-QTL) (2023-08-22, SHA: `4a13d9386f70583dc2bde15477cd313659f48eba`) are documented in this file. 

Note that this branch is in beta and version 2.0.0 is not yet ready for release.


## [2.0.0] - eQTL mapping - 2023-09-21

#### Additions
- Dynamic time usage per rule
- Split `derive_expression_matrices` into `create_expression_matrices` and `prepare_eqtl_input` in order to be able to process all sampels, even those missing genotype information
- Add settings to turn on / off certain calculations
- Added merging of the pools and related rule `barcode_qc_plots` from WG1
- Add optional rule `minimal_postprocess` to output all summary statistics

#### Fixes
- 

#### Changes
- Merged `Create_wp3_inputs` and `Qtl_Snakefile` into one snakemake
- Refactor code to (mostly) PEP8
- Improved clarity of rules by making input / output of rules more explicit
- Moved scripts outside of image for flexibility
- Moved all filepaths to settings file
- Moved all settings to settings file
- Made `ancestry` and `cell level` dynamic
- Refactor to enable parallel preprocessing of pools
- Added `create_seurat` and refactor code to create seurat objects with all metadata files from WG1 and WG2 for memory efficiency
- Made `kinship` optional as it is now moved to WG1
- Refactor code to be more explicit on the LD / eQTL chunks it expects
- Split `derived_feature_annotation_and_chunks` into `feature_annotation` and `eqtl_chunks` in order to split eQTL from LD rules
- Moved `QC_functions.R` to this repo
- Refactor to redo PFlog1pPF as well as recalculate expression PCs after sample removal
- Added example `wg2_pairing.csv` and `qc_mad_final.tab` files to this repo