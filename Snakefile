#!/usr/bin/env python
import re
import os

CHROMOSOMES = [str(chr) for chr in range(1, 23)]
# CHROMOSOMES = ["1"]
CELLTYPES = ["L1/AST", "L1/END", "L1/EX", "L1/IN", "L1/MIC", "L1/OLI", "L1/OPC", "L1/PER"]
# CELLTYPES = ["L1/EX"]

# Add trailing /.
if not config["inputs"]["wg1_folder"].endswith("/"):
    config["inputs"]["wg1_folder"] += "/"
if not config["inputs"]["wg2_folder"].endswith("/"):
    config["inputs"]["wg2_folder"] += "/"
if not config["inputs_extra"]["relative_ancestry_split_path"].endswith("/"):
    config["inputs_extra"]["relative_ancestry_split_path"] += "/"
if not config["outputs"]["output_dir"].endswith("/"):
    config["outputs"]["output_dir"] += "/"

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
# rule run_qtl_mapping.
input_files.append(config["outputs"]["output_dir"] + "input/LimixAnnotationFile.txt")
input_files.extend(expand(config["outputs"]["output_dir"] + "input/{ct}.qtlInput.Pcs.txt",ct=CELLTYPES))
input_files.extend(expand(config["outputs"]["output_dir"] + "input/{ct}.qtlInput.txt", ct=CELLTYPES))
input_files.append(config["outputs"]["output_dir"] + "input/sample.kinship")
input_files.append(config["outputs"]["output_dir"] + "input/L1/smf.txt")

# rule make_temporary_files.
input_files.extend(expand(config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen", chr=CHROMOSOMES))
input_files.extend(expand(config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.bgi", chr=CHROMOSOMES))
input_files.extend(expand(config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.z", chr=CHROMOSOMES))
input_files.extend(expand(config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.sample", chr=CHROMOSOMES))
input_files.extend(expand(config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen_master.txt", chr=CHROMOSOMES))

# rule LD_per_window.
input_files.append(expand(config["outputs"]["output_dir"] + "input/WindowsFilesChecks/chr_{chr}_windows_defined.txt", chr=CHROMOSOMES))

# TODO: According to the settings I should already be able to determine which files will be created
#  allowing me to run the whole pipeline in a single run. However this would require me creating the
#  LimixAnnotationFile, ChunkingFile, and LD windows files in the preflight which can take a bit long..
#  need to think about this.
CHUNK_CHR = []
CHUNK_START = []
CHUNK_END = []
if os.path.exists(config["outputs"]["output_dir"] + "input/ChunkingFile.txt"):
    with open(config["outputs"]["output_dir"] + "input/ChunkingFile.txt") as fp:
        for line in fp:
            re_match = re.match(r"([A-Za-z0-9]+):(\d+)-(\d+)",line.strip())
            CHUNK_CHR.append(re_match[1])
            CHUNK_START.append(re_match[2])
            CHUNK_END.append(re_match[3])

    input_files.extend(expand(config["outputs"]["output_dir"] + "{ct}/top_qtl_results_all.txt.gz", ct=CELLTYPES))
    input_files.extend(expand(config["outputs"]["output_dir"] + "{ct}/qtl.h5.tgz", ct=CELLTYPES))
    input_files.extend(expand(config["outputs"]["output_dir"] + "{ct}/qtl.annotation.tgz", ct=CELLTYPES))
    input_files.append(expand(config["outputs"]["output_dir"] + "{ct}/qtl.permutations.tgz", ct=CELLTYPES))
    input_files.append(config["outputs"]["output_dir"] + "LDMatrices.tgz")


rule all:
    input:
        input_files
    output:
        touch(expand(config["outputs"]["output_dir"] + "{ct}/done.txt", ct=CELLTYPES))

# Import individual rules
include: "includes/prepare_input.smk"
include: "includes/run_qtl.smk"