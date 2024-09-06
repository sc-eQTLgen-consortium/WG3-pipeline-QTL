#!/usr/bin/env python
import re


rule eqtl_chunks:
    input:
        feature_file = config["outputs"]["output_dir"] + "annot_input/{ancestry}/LimixAnnotationFile.txt"
    output:
        chunks = config["outputs"]["output_dir"] + "annot_input/{ancestry}/eQTLChunkingFile.txt"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["eqtl_chunks_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["eqtl_chunks_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["run_qtl"]["eqtl_chunks_time"]]
    threads: config["run_qtl"]["eqtl_chunks_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/derived_eqtl_chunks.R",
        n_genes = config["run_qtl_extra"]["eqtl_chunks_n_genes"],
        out = config["outputs"]["output_dir"] + "annot_input/{ancestry}/"
    log: config["outputs"]["output_dir"] + "log/eqtl_chunks.{ancestry}.log"
    shell:
        """
         singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --feature_file {input.feature_file} \
            --n_genes {params.n_genes} \
            --out_dir {params.out} > {log} 2>&1
        """


rule qtl_kinship_pca:
    input:
        kinship = get_kinship_file,
        wg1_psam = config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_psam"],
        smf = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.smf.txt"
    output:
        pca = config["outputs"]["output_dir"] + "kinship/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.kinship.Pcs.txt",
        pca_var = config["outputs"]["output_dir"] + "kinship/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.kinship.Pcs.var.txt",
        pca_plot = report(config["outputs"]["output_dir"] + "kinship/{ancestry}/QC_figures/{cell_level}.{cell_type}.{qc}.kinship.Pcs.png", category="kinship QC", subcategory="{ancestry} - {cell_level} - {cell_type} - {qc}", caption=config["inputs"]["repo_dir"] + "report_captions/kinship_pca.rst"),
        pca_var_plot = report(config["outputs"]["output_dir"] + "kinship/{ancestry}/QC_figures/{cell_level}.{cell_type}.{qc}.kinship.Pcs.var.png", category="kinship QC", subcategory="{ancestry} - {cell_level} - {cell_type} - {qc}", caption=config["inputs"]["repo_dir"] + "report_captions/kinship_pca_var.rst"),
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["kinship_pca_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["kinship_pca_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["prepare_input"]["kinship_pca_time"]]
    threads: config["prepare_input"]["kinship_pca_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/kinship_pca.R",
        data_out = config["outputs"]["output_dir"] + "kinship/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.kinship.",
        plot_out = config["outputs"]["output_dir"] + "kinship/{ancestry}/QC_figures/{cell_level}.{cell_type}.{qc}.kinship."
    log: config["outputs"]["output_dir"] + "log/kinship_pca.{ancestry}.{cell_level}.{cell_type}.{qc}.log"
    shell:
        """        
         singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --kinship {input.kinship} \
            --wg1_psam {input.wg1_psam} \
            --smf {input.smf} \
            --data_out {params.data_out} \
            --plot_out {params.plot_out} > {log} 2>&1
        """



# This functions takes a normalised expression Seurat object and applies to eQTL mapping filters and prepres the input files required for limix eQTL mapping.
rule prepare_eqtl_input:
    input:
        seurat = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.Qced.Normalized.SCs.Rds"
    output:
        smf = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.smf.txt",
        exp = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.Exp.txt",
        covariates = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.covariates.txt",
        pf = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.qtlInput.txt",
        cf = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.qtlInput.Pcs.txt",
        var = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.qtlInput.Pcs.var.txt",
        pca_plot = report(config["outputs"]["output_dir"] + "expression_input/{ancestry}/QC_figures/{cell_level}.{cell_type}.{qc}.qtlInput.Pcs.png", category="expression QC", subcategory="{ancestry} - {cell_level} - {cell_type} - {qc}", caption=config["inputs"]["repo_dir"] + "report_captions/expression_pca.rst"),
        pca_var_plot = report(config["outputs"]["output_dir"] + "expression_input/{ancestry}/QC_figures/{cell_level}.{cell_type}.{qc}.qtlInput.Pcs.var.png", category="expressiojn QC", subcategory="{ancestry} - {cell_level} - {cell_type} - {qc}", caption=config["inputs"]["repo_dir"] + "report_captions/expression_pca_var.rst"),
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["prepare_eqtl_input_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["prepare_eqtl_input_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["run_qtl"]["prepare_eqtl_input_time"]]
    threads: config["run_qtl"]["prepare_eqtl_input_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/prepare_eqtl_input.R",
        aggregate_fun = config["quality_control_extra"]["prepare_eqtl_input_aggregate_fun"],
        min_n_cells = config["quality_control_extra"]["prepare_eqtl_input_min_n_cells"],
        min_n_samples = config["quality_control_extra"]["prepare_eqtl_input_min_n_samples"],
        individual_aggregate = config["settings_extra"]["individual_aggregate"],
        sample_aggregate = config["settings_extra"]["sample_aggregate"],
        data_out = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.",
        plot_out = config["outputs"]["output_dir"] + "expression_input/{ancestry}/QC_figures/{cell_level}.{cell_type}.{qc}.",
    log: config["outputs"]["output_dir"] + "log/prepare_eqtl_input.{ancestry}.{cell_level}.{qc}.{cell_type}.log"
    shell:
        """
         singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --seurat {input.seurat} \
            --aggregate_fun {params.aggregate_fun} \
            --min_n_cells {params.min_n_cells} \
            --min_n_samples {params.min_n_samples} \
            --individual_aggregate {params.individual_aggregate} \
            --sample_aggregate {params.sample_aggregate} \
            --data_out {params.data_out} \
            --plot_out {params.plot_out} > {log} 2>&1
        """


# Generate a sample QC plot looking at cell count, cell fractions, expression PCs, and kinship PCs.
# Samples are filtered based MAD thresholds selected by the user.
# TODO: consider checking the covariates file as well?
rule qtl_sample_qc:
    input:
        metadata = config["outputs"]["output_dir"] + "expression_input/metadata.tsv.gz",
        # covariates  =  config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.covariates.txt",
        principal_components = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.qtlInput.Pcs.txt",
        full_kinship = config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}.kinship.Pcs.txt",
        kinship = config["outputs"]["output_dir"] + "kinship/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.kinship.Pcs.txt",
        smf = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.smf.txt",
    output:
        fig = report(config["outputs"]["output_dir"] + "expression_input/{ancestry}/QC_figures/{cell_level}.{cell_type}.{qc}.sample_qc.png", category="sample QC", subcategory="{ancestry} - {cell_level} - {cell_type} -  -{qc}", caption=config["inputs"]["repo_dir"] + "report_captions/sample_qc.rst"),
        sample_qc = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/manual_selection/{cell_type}_{qc}_sample_qc.txt",
        pass_qc_samples = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/manual_selection/{cell_type}_{qc}_exclude_smf.txt",
        thresholds = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/manual_selection/{cell_type}_{qc}_threshold_selection.txt"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["qtl_sample_qc_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["qtl_sample_qc_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["run_qtl"]["qtl_sample_qc_time"]]
    threads: config["run_qtl"]["qtl_sample_qc_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/qtl_sample_qc.py",
        donor_cell_threshold = "--donor_cell_threshold " + str(config["quality_control_extra"]["qtl_sample_qc_donor_cell_threshold"]) if config["quality_control_extra"]["qtl_sample_qc_donor_cell_threshold"] is not None else "",
        counts_threshold = "--counts_threshold " + str(config["quality_control_extra"]["qtl_sample_qc_counts_threshold"]) if config["quality_control_extra"]["qtl_sample_qc_counts_threshold"] is not None else "",
        kinship_threshold = "--kinship_threshold " + str(config["quality_control_extra"]["qtl_sample_qc_kinship_threshold"]) if config["quality_control_extra"]["qtl_sample_qc_kinship_threshold"] is not None else "",
        individual_aggregate = config["settings_extra"]["individual_aggregate"],
        sample_aggregate = config["settings_extra"]["sample_aggregate"],
        data_out = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/manual_selection/{cell_type}_{qc}_",
        plot_out = config["outputs"]["output_dir"] + "expression_input/{ancestry}/QC_figures/{cell_level}.{cell_type}.{qc}."
    log: config["outputs"]["output_dir"] + "log/qtl_sample_qc.{ancestry}.{cell_level}.{cell_type}.{qc}.log"
    shell:
        """
         singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --metadata {input.metadata} \
            --ancestry {wildcards.ancestry} \
            --cell_level {wildcards.cell_level} \
            --cell_type {wildcards.cell_type} \
            --qc {wildcards.qc} \
            {params.donor_cell_threshold} \
            --principal_components {input.principal_components} \
            {params.counts_threshold} \
            --full_kinship {input.full_kinship} \
            --kinship {input.kinship} \
            {params.kinship_threshold} \
            --smf {input.smf} \
            --individual_aggregate {params.individual_aggregate} \
            --sample_aggregate {params.sample_aggregate} \
            --data_out {params.data_out} \
            --plot_out {params.plot_out} > {log} 2>&1
        """



rule subset_covariates_file:
    input:
        data = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.qtlInput.Pcs.txt"
    output:
        out = temp(config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.qtlInput.{n_pc}Pcs.txt")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["subset_covariates_file_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["subset_covariates_file_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["run_qtl"]["subset_covariates_file_time"]]
    threads: config["run_qtl"]["subset_covariates_file_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/subset_file.py"
    log: config["outputs"]["output_dir"] + "log/subset_covariates_file.{ancestry}.{cell_level}.{cell_type}.{qc}.{n_pc}.log"
    shell:
        """
         singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --data {input.data} \
            --n_columns {wildcards.n_pc} \
            --out {output.out} > {log} 2>&1
        """


def parse_qtl_chunk(wildcards):
    match = re.match(QTL_CHUNK_PATTERN, wildcards.qtl_chunk)
    return match.group(1), match.group(2), match.group(3)


def get_qtl_bgen(wildcards):
    chr, _, _ = parse_qtl_chunk(wildcards)
    return config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.bgen".format(ancestry=wildcards.ancestry, chr=chr)


def get_qtl_bgen_metafile(wildcards):
    return get_qtl_bgen(wildcards=wildcards) + ".metafile"

def get_qtl_genomic_range(wilcards):
    chr, start, end = parse_qtl_chunk(wilcards)
    return "{chr}:{start}-{end}".format(chr=chr, start=start, end=end)


# Note: add -vf if you want to only test certain variants
rule run_qtl_mapping:
    priority: 10
    input:
        chunks = config["outputs"]["output_dir"] + "annot_input/{ancestry}/eQTLChunkingFile.txt",
        bgen = get_qtl_bgen,
        bgen_metafile = get_qtl_bgen_metafile,
        annotation_file = config["outputs"]["output_dir"] + "annot_input/{ancestry}/LimixAnnotationFile.txt",
        phenotype_file = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/PostQC/{cell_type}.qtlInput.txt",
        covariates_file = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/PostQC/{cell_type}.qtlInput." + str(config["run_qtl_extra"]["n_expression_pcs"]) + "Pcs.txt",
        randomeff_files = get_kinship_file,
        sample_mapping_file = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/PostQC/{cell_type}.smf.txt"
    output:
        h5 = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/qtl_results_{qtl_chunk}.h5",
        # permutation = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/Permutation_information_{feature}.pickle.gz",
        # snp_qc_metrics = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/snp_qc_metrics_naContaining_feature_{feature}.txt.gz",
        snp_metadata = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/snp_metadata_{qtl_chunk}.txt.gz",
        feature_metadata = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/feature_metadata_{qtl_chunk}.txt.gz",
        done = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/{qtl_chunk}.done"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["run_qtl_mapping_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["run_qtl_mapping_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["run_qtl"]["run_qtl_mapping_time"]]
    threads: config["run_qtl"]["run_qtl_mapping_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        python = "/opt/conda/envs/py38/bin/python3.8", # Used to be: python
        script = config["inputs"]["repo_dir"] + "Limix_QTL/run_QTL_analysis_metaAnalysis.py",
        bgen = lambda wildcards, input: input.bgen.rstrip(".bgen"),
        output_directory = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/",
        window = config["run_qtl_extra"]["window_size"],
        genomic_range = get_qtl_genomic_range,
        covariates_file = lambda wildcards, input: "--covariates_file " + str(input.covariates_file) if input.covariates_file else "",
        maf = config["run_qtl_extra"]["maf"],
        hwe = config["run_qtl_extra"]["hwe"],
        call_rate = config["run_qtl_extra"]["call_rate"],
        block_size = config["run_qtl_extra"]["block_size"],
        number_of_permutations = config["run_qtl_extra"]["n_permutations"],
        relatedness_score = config["run_qtl_extra"]["relatedness_score"],
        minimum_test_samples = config["run_qtl_extra"]["minimum_test_samples"],
    log: config["outputs"]["output_dir"] + "log/run_qtl_mapping.{ancestry}.{cell_level}.{cell_type}.chr_{qtl_chunk}.log"
    shell:
        """
         singularity exec --bind {params.bind} {params.sif} {params.python} {params.script} \
            --bgen {params.bgen} \
            --annotation_file {input.annotation_file} \
            --phenotype_file {input.phenotype_file} \
            --output_directory {params.output_directory} \
            --window {params.window} \
            --genomic_range {params.genomic_range} \
            {params.covariates_file} \
            --randomeff_files {input.randomeff_files} \
            --sample_mapping_file {input.sample_mapping_file} \
            --minor_allel_frequency {params.maf} \
            --hardy_weinberg_equilibrium {params.hwe} \
            --call_rate {params.call_rate} \
            --block_size {params.block_size} \
            --number_of_permutations {params.number_of_permutations} \
            --relatedness_score {params.relatedness_score} \
            --minimum_test_samples {params.minimum_test_samples} \
            --gaussianize_method gaussnorm \
            --cis > {log} 2>&1
        
         singularity exec --bind {params.bind} {params.sif} touch {output.h5} {output.snp_metadata} {output.feature_metadata} {output.done}
        """


rule minimal_postprocess:
    input:
        # h5 = lambda wildcards: expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/qtl_results_{qtl_chunk}.h5", qtl_chunk=QTL_CHUNKS, allow_missing=True),
        # snp_metadata = lambda wildcards: expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/snp_metadata_{qtl_chunk}.txt.gz", qtl_chunk=QTL_CHUNKS, allow_missing=True),
        # feature_metadata = lambda wildcards: expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/feature_metadata_{qtl_chunk}.txt.gz", qtl_chunk=QTL_CHUNKS, allow_missing=True)
        done = lambda wildcards: expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/{qtl_chunk}.done", qtl_chunk=QTL_CHUNKS, allow_missing=True),
    output:
        out = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/{resultsfile}.txt.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["minimal_postprocess_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["minimal_postprocess_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["run_qtl"]["minimal_postprocess_time"]]
    threads: config["run_qtl"]["minimal_postprocess_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        python = "/opt/conda/envs/py38/bin/python3.8", # Used to be: python
        script = config["inputs"]["repo_dir"] + "Limix_QTL/post_processing/minimal_postprocess.py",
        input_dir = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/",
        ouput_dir = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/",
        minimal_reporting_p = 1.0,
        minimal_reporting_featureP = 1.0,
        top_feature_based = lambda wildcards: "--top_feature_based" if wildcards.resultsfile.startswith("top_") else ""
    log: config["outputs"]["output_dir"] + "log/minimal_postprocess.{ancestry}.{cell_level}.{cell_type}.{resultsfile}.log"
    shell:
        """
         singularity exec --bind {params.bind} {params.sif} {params.python} {params.script} \
            --input_dir {params.input_dir} \
            --ouput_dir {params.ouput_dir} \
            --single_file_output \
            --write_compressed \
            --minimimal_reporting_p {params.minimal_reporting_p} \
            --minimal_reporting_featureP {params.minimal_reporting_featureP} \
            {params.top_feature_based} > {log} 2>&1
        """


rule combine_qtl_mapping_runtime:
    input:
        log = lambda wildcards: expand(config["outputs"]["output_dir"] + "log/run_qtl_mapping.{ancestry}.{cell_level}.{cell_type}.chr_{qtl_chunk}.log", qtl_chunk=QTL_CHUNKS, allow_missing=True)
    output:
        runtime = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/runtime.tsv.gz",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["combine_qtl_mapping_runtime_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["combine_qtl_mapping_runtime_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["run_qtl"]["combine_qtl_mapping_runtime_time"]]
    threads: config["run_qtl"]["combine_qtl_mapping_runtime_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/combine_qtl_mapping_runtime.py",
        input = config["outputs"]["output_dir"] + "log/run_qtl_mapping.{ancestry}.{cell_level}.{cell_type}"
    log: config["outputs"]["output_dir"] + "log/combine_qtl_mapping_runtime.{ancestry}.{cell_level}.{cell_type}.log"
    shell:
        """
         singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --input {params.input} \
            --output {output.runtime} > {log} 2>&1
        """


rule multiple_testing_correction:
    input:
        top_hits = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/top_qtl_results_all.txt.gz"
    output:
        top_hits = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/top_qtl_results_all_multest.txt",
        fig1 = report(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/QC_figures/{cell_type}.top_qtl_results_all_overview.png", category="Multiple testing correction", subcategory="{ancestry} - {cell_level} - {cell_type}", caption=config["inputs"]["repo_dir"] + "report_captions/qvalues_overview.rst"),
        fig2 = report(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/QC_figures/{cell_type}.top_qtl_results_all_hist.png", category="Multiple testing correction", subcategory="{ancestry} - {cell_level} - {cell_type}", caption=config["inputs"]["repo_dir"] + "report_captions/qvalues_hist.rst")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["multiple_testing_correction_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["multiple_testing_correction_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["run_qtl"]["multiple_testing_correction_time"]]
    threads: config["run_qtl"]["multiple_testing_correction_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/multiple_testing_correction.R",
        pvalue = "empirical_feature_p_value",
        qvalue = "empirical_feature_q_value",
        alpha = 0.05,
        data_out = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/top_qtl_results_all",
        plot_out = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/QC_figures/{cell_type}.top_qtl_results_all",
        suffix = "_multest"
    log: config["outputs"]["output_dir"] + "log/multiple_testing_correction.{ancestry}.{cell_level}.{cell_type}.log"
    shell:
        """
         singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --input {input.top_hits} \
            --pvalue {params.pvalue} \
            --qvalue {params.qvalue} \
            --alpha {params.alpha} \
            --data_out {params.data_out} \
            --plot_out {params.plot_out} \
            --suffix {params.suffix} > {log} 2>&1
        """


# TODO: remove 'qtl' folder after compression.
rule compress_qtl:
    input:
        # h5 = expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/qtl_results_{qtl_chunk}.h5", qtl_chunk=QTL_CHUNKS, allow_missing=True),
        # # permutation = expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/Permutation_information_{feature}.pickle.gz", feature=get_features, allow_missing=True),
        # # snp_qc_metrics = expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/snp_qc_metrics_naContaining_feature_{feature}.txt.gz", feature=get_features, allow_missing=True),
        # snp_metadata = expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/snp_metadata_{qtl_chunk}.txt.gz", qtl_chunk=QTL_CHUNKS, allow_missing=True),
        # feature_metadata = expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/feature_metadata_{qtl_chunk}.txt.gz", qtl_chunk=QTL_CHUNKS, allow_missing=True),
        done = lambda wildcards: expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/{qtl_chunk}.done", qtl_chunk=QTL_CHUNKS, allow_missing=True),
    output:
        h5 = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl.h5.tgz",
        anno = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl.annotation.tgz",
        perm = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl.permutations.tgz",
        done = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/compress_qtl.done"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["compress_qtl_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["compress_qtl_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["run_qtl"]["compress_qtl_time"]]
    threads: config["run_qtl"]["compress_qtl_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        qtl_output_dir = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/",
        h5_merge_list = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/h5_merge_list.txt",
        feature_merge_list = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/feature_merge_list.txt",
        anno_merge_list = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/anno_merge_list.txt"
    log: config["outputs"]["output_dir"] + "log/compress_qtl.{ancestry}.{cell_level}.{cell_type}.log"
    shell:
        """
         singularity exec --bind {params.bind} {params.sif} find {params.qtl_output_dir} -type f -name \*.h5 > {params.h5_merge_list}
         singularity exec --bind {params.bind} {params.sif} tar -czf {output.h5} -T {params.h5_merge_list}

         singularity exec --bind {params.bind} {params.sif} find {params.qtl_output_dir} -type f -name \*.pickle.gz > {params.feature_merge_list}
         singularity exec --bind {params.bind} {params.sif} tar -czf {output.perm} -T {params.feature_merge_list}
        
         singularity exec --bind {params.bind} {params.sif} find {params.qtl_output_dir} -type f -name \*.txt.gz > {params.anno_merge_list}
         singularity exec --bind {params.bind} {params.sif} tar -cf {output.anno} -T {params.anno_merge_list}
        
         singularity exec --bind {params.bind} {params.sif} touch {output.done}
        """
