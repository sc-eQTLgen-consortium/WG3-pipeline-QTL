#!/usr/bin/env python
import re
import os

CHROMOSOMES = [str(chr) for chr in range(1, 23)]
CELLTYPES = ["L1/AST", "L1/END", "L1/EX", "L1/IN", "L1/MIC", "L1/OLI", "L1/OPC", "L1/PER"]


# Add trailing /.
if not config["inputs"]["wg1_folder"].endswith("/"):
    config["inputs"]["wg1_folder"] += "/"
if not config["inputs"]["wg2_folder"].endswith("/"):
    config["inputs"]["wg2_folder"] += "/"
if not config["inputs_extra"]["relative_ancestry_split_path"].endswith("/"):
    config["inputs_extra"]["relative_ancestry_split_path"] += "/"

# Check if the input files exist.
if not os.path.exists(config["refs"]["gtf_annotation_file"]):
    logger.info("Error, the GTF annotation file '{}' does not exist.\n\nExiting.".format(config["refs"]["gtf_annotation_file"]))
    exit("MissingGTFAnnotationInput")
if not os.path.exists(config["inputs"]["unimputed_folder"] + ".pgen") or not os.path.exists(config["inputs"]["unimputed_folder"] + ".pvar") or not os.path.exists(config["inputs"]["unimputed_folder"] + ".psam"):
    logger.info("Error, the unimputed genotype file '{}'[.pgen/.pvar/.psam] does not exist.\n\nExiting.".format(config["inputs"]["unimputed_folder"]))
    exit("MissinGenotypeInput")
if not os.path.exists(config["inputs"]["cell_annotation"]):
    logger.info("Error, the cell annotations file '{}' does not exist.\n\nExiting.".format(config["inputs"]["cell_annotation"]))
    exit("MissingCellAnnotationInput")
for relative_path in ["relative_wg1_psam", "relative_wg1_singlets_assigned"]:
    if not os.path.exists(config["inputs"]["wg1_folder"] + config["inputs_extra"][relative_path]):
        logger.info("Could not find the {} file. Please check that the file exists.\n\nExiting.".format(config["inputs"]["wg1_folder"] + config["inputs_extra"][relative_path]))
        exit("MissingWG1InputFile")

#####################
######## ALL ########
#####################
input_files = []

# Add the endpoints for pre_processing.
input_files.extend(expand(config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.metafile",chr=CHROMOSOMES))
input_files.append(config["outputs"]["output_dir"] + "input/AllMetaData.debug.txt")
input_files.append(config["outputs"]["output_dir"] + "input/sample.kinship")
input_files.append(config["outputs"]["output_dir"] + "input/ChunkingFile.txt")
input_files.append(config["outputs"]["output_dir"] + "input/LimixAnnotationFile.txt")
input_files.append(config["outputs"]["output_dir"] + "input/smf.txt")
input_files.extend(config["outputs"]["output_dir"] + "input/WindowsFilesChecks/chr_{chr}_windows_defined.txt", chr=CHROMOSOMES)

# TODO: According to the settings I should already be able to determine which files will be created
#  allowing me to run the whole pipeline in a single run.
if os.path.exists(config["inputs"]["wg3_folder"] + "input/ChunkingFile.txt"):
    chunk_chr, chunk_start, chunk_end = [], [], []
    with open(config["inputs"]["wg3_folder"] + "input/ChunkingFile.txt") as fp:
        for line in fp:
            re_match = re.match(r"([A-Za-z0-9]+):(\d+)-(\d+)",line.strip())
            chunk_chr.append(re_match[1])
            chunk_start.append(re_match[2])
            chunk_end.append(re_match[3])

    input_files.extend(expand(config["outputs"]["output_dir"] + "{ct}/top_qtl_results_all.txt.gz", ct=CELLTYPES))
    input_files.extend(expand(config["outputs"]["output_dir"] + "{ct}/qtl.h5.tgz", ct=CELLTYPES))
    input_files.extend(expand(config["outputs"]["output_dir"] + "{ct}/qtl.annotation.tgz", ct=CELLTYPES))
    input_files.append(expand(config["outputs"]["output_dir"] + "{ct}/qtl.permutations.tgz", ct=CELLTYPES))
    input_filesa.append(config["outputs"]["output_dir"] + "LDMatrices.tgz")


rule all:
    input:
        input_files
    output:
        touch(expand(config["outputs"]["output_dir"] + "{ct}/done.txt", ct=CELLTYPES))

# Import individual rules
include: "includes/prepare_input.smk"
include: "includes/run_qtl.smk"