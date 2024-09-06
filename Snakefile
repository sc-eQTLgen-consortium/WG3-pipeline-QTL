#!/usr/bin/env python
import pandas as pd
import gzip
import re
import os

CHROMOSOMES = [str(chr) for chr in range(1, 23)]

CHR_PATTERN = "[0-9]{1,2}|X|Y|MT"
QTL_CHUNK_PATTERN = "([0-9]{1,2}|X|Y|MT)_([0-9]+)_([0-9]+)"
LD_CHUNK_PATTERN = "chr_([0-9]{1,2}|X|Y|MT)_window_([0-9]+)_([0-9]+)_([0-9]+)"

# Add trailing /.
if not config["inputs"]["repo_dir"].endswith("/"):
    config["inputs"]["repo_dir"] += "/"
if not config["inputs"]["wg1_genotype_folder"].endswith("/"):
    config["inputs"]["wg1_genotype_folder"] += "/"
if not config["inputs"]["wg1_demultiplex_folder"].endswith("/"):
    config["inputs"]["wg1_demultiplex_folder"] += "/"
if not config["inputs"]["wg2_folder"].endswith("/"):
    config["inputs"]["wg2_folder"] += "/"
if not config["refs"]["ref_dir"].endswith("/"):
        config["refs"]["ref_dir"] += "/"
if not config["outputs"]["output_dir"].endswith("/"):
    config["outputs"]["output_dir"] += "/"

# Function that quits the program.
# This is seperated so I can ignore errors if I want to.
if config["settings_extra"]["ignore_file_checks"]:
    logger.warning("Ignoring errors, please note that this may cause rules to crash!\n")

def stop(message):
    if config["settings_extra"]["ignore_file_checks"]:
        return
    exit(message)

# Check if the input files exist.
if not os.path.exists(config["inputs"]["singularity_image"]):
    logger.info("Error, the singularity image does not exist.\n\nExiting.")
    stop("MissingSIFFile")
if config["settings"]["calculate_qtl"] and not os.path.exists(config["inputs"]["poolsheet_path"]):
    logger.info("Error, the poolsheet file does not exist.\n\nExiting.")
    stop("MissingPoolSheetFile")
if (config["settings"]["calculate_qtl"] or config["settings"]["calculate_ld"]) and not os.path.exists(config["refs"]["gtf_annotation_file"]):
    logger.info("Error, the GTF annotation file '{}' does not exist.\n\nExiting.".format(config["refs"]["gtf_annotation_file"]))
    stop("MissingGTFAnnotationInput")

genotype_vcf = (config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_imputed_genotype_vcf"]).replace("{ancestry}", config["settings"]["ancestry"])
if (config["settings"]["calculate_qtl"] or config["settings"]["calculate_ld"]) and not os.path.exists(genotype_vcf):
    logger.info("Error, the genotype VCF file '{}' does not exist.\n\nExiting.".format(genotype_vcf))
    stop("MissingGenotypeInput")

psam_file = config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_psam"]
if not os.path.exists(psam_file):
    logger.info("Error, the PSAM file '{}' does not exist.\n\nExiting.".format(psam_file))
    stop("MissingPSAMInput")

het_filtered_vcf = (config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_het_filtered"]).replace("{ancestry}", config["settings"]["ancestry"])
kinship_file = (config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_kinship"]).replace("{ancestry}", config["settings"]["ancestry"])
if config["settings"]["calculate_qtl"] and not os.path.exists(het_filtered_vcf) and not os.path.exists(kinship_file):
    logger.info("Could not find the kinship ('{}') nor the HET filtered VCF ('{}') file. Please check that one of the files exists.\n\nExiting.".format(het_filtered_vcf, kinship_file))
    stop("MissingKinshipInputFile")

if not os.path.exists(config["inputs"]["wg1_demultiplex_folder"] + config["inputs_extra"]["relative_wg1_metadata"]):
    logger.info("Could not find the {} file. Please check that the file exists.\n\nExiting.".format(config["inputs"]["wg1_demultiplex_folder"] + config["inputs_extra"]["relative_wg1_metadata"]))
    stop("MissingWG1MetadataFile")
if not os.path.exists(config["inputs"]["wg2_folder"] + config["inputs_extra"]["relative_wg2_metadata"]):
    logger.info("Could not find the {} file. Please check that the file exists.\n\nExiting.".format(config["inputs"]["wg2_folder"] + config["inputs_extra"]["relative_wg2_metadata"]))
    stop("MissingWG2MetadataFile")
if not os.path.exists(config["inputs"]["cell_annotation_file"]):
    logger.info("Error, the cell annotation file '{}' does not exist.\n\nExiting.".format(config["inputs"]["cell_annotation_file"]))
    stop("MissingCellAnnotationFile")

# Check the rb / mt genes for plotting.
for genes_file in ["ribosomal_genes", "mitochondrial_genes", "qc_mad"]:
    if not os.path.exists(config["refs"]["ref_dir"] + config["refs_extra"][genes_file]):
        logger.info("Could not find the {} file. Please check that the file exists.\n\nExiting.".format(config["refs"]["ref_dir"] + config["refs_extra"][genes_file]))
        stop("MissingReferenceFile")


def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)


logger.info("Loading the input poolsheet")
POOL_DF = pd.read_csv(config["inputs"]["poolsheet_path"], sep="\t", dtype=str)
POOL_DF.fillna("NA", inplace=True)
POOL_DF.index = POOL_DF["Pool"]

if "Pool" not in POOL_DF.columns:
    logger.info("\tError, missing 'Pool' column in poolsheet file for the selected methods.\n\nExiting.")
    stop("InvalidPoolSheetFile")
if not POOL_DF["Pool"].is_unique:
    logger.info("\tError, your 'Pool' column contains duplicates, please make sure all values are unique.\n\nExiting.")
    stop("InvalidPoolSheetFile")

logger.info("\tValid.")
POOLS = POOL_DF["Pool"].tolist()

required_columns = {
    'Barcode': (True, False, {}),
    'sequencing_platform': (False, False, {}),
    'sequencing_run': (False, False, {}),
    'sequencing_lane': (False, False, {}),
    'scrna_platform': (False, False, {}),
    'plate_based': (False, False, {"Y": "Yes", "N": "No", "NONE": "unknown"}),
    'umi_based': (False, False, {"Y": "Yes", "N": "No", "NONE": "unknown"}),
    'biomaterial': (False, False, {}),
    'sorting': (False, False, {}),
    'cell_treatment': (False, False, {}),
    'sample_condition': (False, False, {})
}

# Loading the cell annotation.
logger.info("Validating cell annotations:")
cell_annot_df = pd.read_csv(config["inputs"]["cell_annotation_file"], sep="\t", dtype=str, keep_default_na=False, nrows=1)

# Check for missing columns.
missing_columns = [column for column in required_columns.keys() if not column in cell_annot_df.columns]
if len(missing_columns) > 0:
    logger.info("\tThe column names of your cell annotation file are not correct.\n\t"
                "The columns that you are missing or whose spelling does not match the required input is/are: {}.\n\t"
                "They should be: 'Barcode', 'sequencing_platform', 'sequencing_run', 'sequencing_lane', 'scrna_platform', 'plate_based','umi_based', 'biomaterial', 'sorting','cell_treatment', 'sample_condition'.\n\t"
                "If the names look the same, check that the file is tab separated, without any spaces or other weird characters.".format(",".join(missing_columns)))
    stop("InvalidCellAnnotation")

# Check if the contents is valid.
cell_annot_is_valid = True
for column, (must_be_unique, na_allowed, valid_values) in required_columns.items():
    if must_be_unique and not cell_annot_df[column].is_unique:
        logger.info("\tYour {} column contains duplicates, please make sure all values are unique.".format(column))
        cell_annot_is_valid = False

    # Check if there are NONE in the column.
    if not na_allowed and "NONE" in cell_annot_df[column]:
        logger.info("\tYour {} column is missing entries. NONE values are not allowed in this column.\n"
                    "Please make sure that all contents of the study column have a string entry.".format(column))
        cell_annot_is_valid = False

    # Check if the values correspond with the accepted input values.
    if valid_values and not cell_annot_df[column].isin(valid_values.keys()).all():
        logger.info("\tYour {} column does not have just {}, please make sure all values in this column are {}.".format(column,", ".join([str(x) for x in valid_values]),", ".join(["{} ({})".format(key,value) for key, value in valid_values.items()])))
        cell_annot_is_valid = False

if not cell_annot_is_valid:
    logger.info("\n\nExiting.")
    stop("InvalidCellAnnotation")

logger.info("\tValid.")

#####################
######## ALL ########
#####################
input_files = []

# Create the expression matrices for gene networks.
if config["settings"]["save_all_samples"]:
    logger.info("Combining all RNA-seq samples for downstream analyses.")
    for cell_type in config["settings"]["cell_types"]:
        input_files.append(config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.Qced.Normalized.SCs.Rds".format(
            ancestry="ALL",
            qc="PreQC",
            cell_level=config["settings"]["cell_level"],
            cell_type=cell_type))

QTL_CHUNKS = []
if config["settings"]["calculate_qtl"]:
    logger.info("Calculating eQTLs.")

    # First create the chunks.
    eqtl_chunks_path = config["outputs"]["output_dir"] + "annot_input/{ancestry}/eQTLChunkingFile.txt".format(ancestry=config["settings"]["ancestry"])
    if os.path.exists(eqtl_chunks_path):
        with gzopen(eqtl_chunks_path,mode='r') as f:
            for i, line in enumerate(f):
                match = re.match("([0-9]{1,2}|X|Y|MT):([0-9]+)-([0-9]+)",line)
                chunk = "{chr}_{start}_{end}".format(chr=match.group(1),start=match.group(2),end=match.group(3))
                QTL_CHUNKS.append(chunk)
        f.close()
    else:
        logger.info("Creating eQTL chunks.")
        input_files.append(eqtl_chunks_path)

    # Create the barcode QC figures.
    input_files.append(config["outputs"]["output_dir"] + "expression_input/QC_figures/qc_plots.done")

    if len(QTL_CHUNKS) == 0:
        print("\tError, no eQTL chunks")
        stop("InvalideQTLChunks")

    for cell_type in config["settings"]["cell_types"]:
        logger.info("  Cell level: {}, cell type: {}".format(config["settings"]["cell_level"],cell_type))

        # Get the sample QC threshold selections.
        threshold_select_path = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/manual_selection/{cell_type}_{qc}_threshold_selection.txt".format(
            ancestry=config["settings"]["ancestry"],
            cell_level=config["settings"]["cell_level"],
            cell_type=cell_type,
            qc="PreQC")
        input_files.append(threshold_select_path)

        # Parse the sample QC threshold selection.
        if os.path.exists(threshold_select_path):
            threshold_select_df = pd.read_csv(threshold_select_path,sep="\t",header=0,index_col=None)
            if not threshold_select_df["PASS"].all():
                logger.info("    Warning, waiting for user to validate QC thresholds. Please make sure all values in '{}' column 'PASS' are True. If you want to select different MAD thresholds, adjust the config file and remove the this file.".format(os.path.basename(threshold_select_path)))
                continue
            for threshold in ["donor_cell_threshold", "counts_threshold", "kinship_threshold"]:
                if str(threshold_select_df.loc[threshold_select_df["Threshold"] == threshold, "Value"].values[0]) != str(config["quality_control_extra"]["qtl_sample_qc_" + threshold]):
                    logger.info("    Error, the threshold for '{}' in the settings file does not match the value in '{}'. If you want to select different MAD thresholds, adjust the config file and remove the this file.\n\nExiting.".format(threshold,os.path.basename(threshold_select_path)))
                    exit("InvalidSampleQCThreshold")

            logger.info("    Sample QC passed.")

            # Create the pass QC output. This forces to rerun the sample QC rule eventhough we do not need that output.
            input_files.append(config["outputs"]["output_dir"] + "expression_input/{ancestry}/QC_figures/{cell_level}.{cell_type}.{qc}.sample_qc.png".format(
                ancestry=config["settings"]["ancestry"],
                cell_level=config["settings"]["cell_level"],
                cell_type=cell_type,
                qc="PostQC"))

            # Create the eQTL output.
            input_files.append(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/top_qtl_results_all_multest.txt".format(
                ancestry=config["settings"]["ancestry"],
                cell_level=config["settings"]["cell_level"],
                cell_type=cell_type))

            # Also combine the runtime.
            input_files.append(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/runtime.tsv.gz".format(
                ancestry=config["settings"]["ancestry"],
                cell_level=config["settings"]["cell_level"],
                cell_type=cell_type))

            # Add the optional eQTL outputs.
            if config["settings"]["output_flat_qtl_results"]:
                input_files.append(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl_results_all.txt.gz".format(
                    ancestry=config["settings"]["ancestry"],
                    cell_level=config["settings"]["cell_level"],
                    cell_type=cell_type))
            if config["settings"]["compress_qtl"]:
                input_files.append(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/compress_qtl.done".format(
                    ancestry=config["settings"]["ancestry"],
                    cell_level=config["settings"]["cell_level"],
                    cell_type=cell_type))
        else:
            logger.info("    Waiting on sample QC results.")

LD_CHUNKS = []
if config["settings"]["calculate_ld"]:
    logger.info("Creating compressed LD files.")

    ld_chunks_path = config["outputs"]["output_dir"] + "annot_input/{ancestry}/LDChunkingFile.txt".format(ancestry=config["settings"]["ancestry"])
    if os.path.exists(ld_chunks_path):
        with gzopen(ld_chunks_path,mode='r') as f:
            for line in f:
                match = re.match("([0-9]{1,2}|X|Y|MT):([0-9]+):([0-9]+)-([0-9]+)",line)
                LD_CHUNKS.append("chr_{chr}_window_{num}_{start}_{end}".format(chr=match.group(1),num=match.group(2),start=match.group(3),end=match.group(4)))
        f.close()
    else:
        logger.info("Creating LD chunks.")
        input_files.append(ld_chunks_path)

    # Create the LD output.
    if len(LD_CHUNKS) == 0:
        print("\tError, no LD chunks")
        stop("InvalidLDChunks")

    if config["settings"]["compress_ld"]:
        input_files.append(config["outputs"]["output_dir"] + "output/{ancestry}/LDMatrices/LDMatrices.tgz".format(
            ancestry=config["settings"]["ancestry"]))
    else:
        for ld_chunk in LD_CHUNKS:
            ld_chr = re.match(LD_CHUNK_PATTERN,ld_chunk).group(1)
            input_files.extend([
                config["outputs"]["output_dir"] + "output/{ancestry}/LDMatrices/chr{ld_chr}/{ld_chunk}_low_dim.pkl.gz".format(ancestry=config["settings"]["ancestry"],ld_chr=ld_chr,ld_chunk=ld_chunk),
                config["outputs"]["output_dir"] + "output/{ancestry}/LDMatrices/chr{ld_chr}/{ld_chunk}_components.pkl.gz".format(ancestry=config["settings"]["ancestry"],ld_chr=ld_chr,ld_chunk=ld_chunk),
                config["outputs"]["output_dir"] + "output/{ancestry}/LDMatrices/chr{ld_chr}/{ld_chunk}_means.pkl.gz".format(ancestry=config["settings"]["ancestry"],ld_chr=ld_chr,ld_chunk=ld_chunk),
                config["outputs"]["output_dir"] + "output/{ancestry}/LDMatrices/chr{ld_chr}/{ld_chunk}.pkl.gz".format(ancestry=config["settings"]["ancestry"],ld_chr=ld_chr,ld_chunk=ld_chunk)
            ])


rule all:
    input:
        input_files


wildcard_constraints:
    ancestry = "\w+",
    pool = "[\w-]+",
    cell_level = "\w+",
    cell_type = "[A-Za-z.]+",
    qc = "\w+",
    chr = CHR_PATTERN,
    qtl_chunk = QTL_CHUNK_PATTERN,
    top = "\w+",
    norm_file = "[A-Za-z.]+",
    ld_chr = CHR_PATTERN,
    ld_chunk = LD_CHUNK_PATTERN


# Import individual rules
include: "includes/prepare_input.smk"
include: "includes/run_qtl.smk"
include: "includes/calculate_ld.smk"
