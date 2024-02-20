#!/usr/bin/env python
import os


# In this rule, we want to read this meta data file and calculate MAD stuff. This has to be over all pools together!
# A tag is added to the metadata to mark if the barcode passed QC or not.
# TODO: parameter md_vars and downsampling is not implemented
rule barcode_qc:
    input:
        poolsheet = config["inputs"]["poolsheet_path"],
        metadata = expand(config["outputs"]["output_dir"] + "expression_input/pools/{pool}.metadata.tsv.gz", pool=POOLS)
    output:
        metadata = config["outputs"]["output_dir"] + "expression_input/metadata.tsv.gz",
        tagged = config["outputs"]["output_dir"] + "expression_input/metadata.tagged.tsv.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["quality_control"]["barcode_qc_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["quality_control"]["barcode_qc_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["quality_control"]["barcode_qc_time"]]
    threads: config["quality_control"]["barcode_qc_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/barcode_qc.R",
        functions_fn = config["inputs"]["repo_dir"] + "scripts/QC_functions.R",
        input_dir = config["outputs"]["output_dir"] + "expression_input/pools/",
        level = config["quality_control_extra"]["barcode_qc_level"],
        qc_mad = config["refs"]["ref_dir"] + config["refs_extra"]["qc_mad"],
        out_dir = config["outputs"]["output_dir"] + "expression_input/"
    log: config["outputs"]["output_dir"] + "log/barcode_qc.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --functions_fn {params.functions_fn} \
            --poolsheet {input.poolsheet} \
            --input_dir {params.input_dir} \
            --level {params.level} \
            --qc_mad {params.qc_mad} \
            --out_dir {params.out_dir}
        """

rule barcode_qc_plots:
    input:
        metadata = config["outputs"]["output_dir"] + "expression_input/metadata.tsv.gz"
    output:
        fig1 = report(config["outputs"]["output_dir"] + "expression_input/QC_figures/nCount_RNA_violin_MAD_All.png", category="expression QC", caption=config["inputs"]["repo_dir"] + "report_captions/QC_plots_nCount.rst"),
        fig2 = report(config["outputs"]["output_dir"] + "expression_input/QC_figures/nCount_RNA_violin_MADper_Pool.png", category="expression QC", caption=config["inputs"]["repo_dir"] + "report_captions/QC_plots_nCount_RNA_MADall.rst"),
        fig3 = report(config["outputs"]["output_dir"] + "expression_input/QC_figures/nCount_RNA_violin_noMADlines.png", category="expression QC", caption=config["inputs"]["repo_dir"] + "report_captions/QC_plots_nCount.rst"),
        fig4 = report(config["outputs"]["output_dir"] + "expression_input/QC_figures/nFeature_RNA_violin_MAD_All.png", category="expression QC", caption=config["inputs"]["repo_dir"] + "report_captions/QC_plots_feature_MADall.rst"),
        fig5 = report(config["outputs"]["output_dir"] + "expression_input/QC_figures/nFeature_RNA_violin_MADper_Pool.png", category="expression QC", caption=config["inputs"]["repo_dir"] + "report_captions/QC_plots_feature_MADperPool.rst"),
        fig6 = report(config["outputs"]["output_dir"] + "expression_input/QC_figures/nFeature_RNA_violin_noMADlines.png", category="expression QC", caption=config["inputs"]["repo_dir"] + "report_captions/QC_plots_feature.rst"),
        fig7 = report(config["outputs"]["output_dir"] + "expression_input/QC_figures/nFeatures_vs_percentMT_QC_scatter_colorPool.png", category="expression QC", caption=config["inputs"]["repo_dir"] + "report_captions/QC_plots_features_mt_pool.rst"),
        fig8 = report(config["outputs"]["output_dir"] + "expression_input/QC_figures/nFeatures_vs_percentMT_QC_scatter_w_MADlines.png", category="expression QC", caption=config["inputs"]["repo_dir"] + "report_captions/QC_plots_features_mt_MAD.rst"),
        fig9 = report(config["outputs"]["output_dir"] + "expression_input/QC_figures/percent.mt_violin_MAD_All.png", category="expression QC", caption=config["inputs"]["repo_dir"] + "report_captions/QC_plots_mt_MADall.rst"),
        fig10 = report(config["outputs"]["output_dir"] + "expression_input/QC_figures/percent.mt_violin_MADper_Pool.png", category="expression QC", caption=config["inputs"]["repo_dir"] + "report_captions/QC_plots_mt_MADpool.rst"),
        fig11 = report(config["outputs"]["output_dir"] + "expression_input/QC_figures/percent.mt_violin_noMADlines.png", category="expression QC", caption=config["inputs"]["repo_dir"] + "report_captions/QC_plots_mt.rst"),
        fig12 = report(config["outputs"]["output_dir"] + "expression_input/QC_figures/percent.rb_violin_MAD_All.png", category="expression QC", caption=config["inputs"]["repo_dir"] + "report_captions/QC_plots_rb_MADall.rst"),
        fig13 = report(config["outputs"]["output_dir"] + "expression_input/QC_figures/percent.rb_violin_MADper_Pool.png", category="expression QC", caption=config["inputs"]["repo_dir"] + "report_captions/QC_plots_rb_MADperPool.rst"),
        fig14 = report(config["outputs"]["output_dir"] + "expression_input/QC_figures/percent.rb_violin_noMADlines.png", category="expression QC", caption=config["inputs"]["repo_dir"] + "report_captions/QC_plots_rb.rst"),
        fig15 = report(config["outputs"]["output_dir"] + "expression_input/QC_figures/UMI_vs_Genes_QC_scatter.png", category="expression QC", caption=config["inputs"]["repo_dir"] + "report_captions/QC_plots_UMI_features.rst"),
        fig16 = report(config["outputs"]["output_dir"] + "expression_input/QC_figures/UMI_vs_Genes_QC_scatter_w_MADlines.png", category="expression QC", caption=config["inputs"]["repo_dir"] + "report_captions/QC_plots_UMI_features_MADall.rst"),
        fig17 = report(config["outputs"]["output_dir"] + "expression_input/QC_figures/UMI_vs_percentMT_QC_scatter_colorPool.png", category="expression QC", caption=config["inputs"]["repo_dir"] + "report_captions/QC_plots_UMI_mt_pool.rst"),
        fig18 = report(config["outputs"]["output_dir"] + "expression_input/QC_figures/UMI_vs_percentMT_QC_scatter_w_MADlines.png", category="expression QC", caption=config["inputs"]["repo_dir"] + "report_captions/QC_plots_UMI_mt_MDA_all.rst"),
        done = config["outputs"]["output_dir"] + "expression_input/QC_figures/qc_plots.done"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["quality_control"]["barcode_qc_plots_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["quality_control"]["barcode_qc_plots_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["quality_control"]["barcode_qc_plots_time"]]
    threads: config["quality_control"]["barcode_qc_plots_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/barcode_qc_plots.R",
        out = config["outputs"]["output_dir"] + "expression_input/QC_figures/"
    log: config["outputs"]["output_dir"] + "log/barcode_qc_plots.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --metadata {input.metadata} \
            --out {params.out}
        singularity exec --bind {params.bind} {params.sif} touch {output.done}
        """


# For each cell type load the pool seurat object, filter on singlet, barcode QC tag, cell type, ancestry (optional), treatment and combine into one seurat object.
# Then apply PF - log - PF: calculate total number of UMI over all samples (meanSampleSum).
# Optional: sample are removed baed on rule 'sample_qc'.
rule create_expression_matrices:
    input:
        poolsheet = config["inputs"]["poolsheet_path"],
        seurat = expand(config["outputs"]["output_dir"] + "expression_input/pools/{pool}.rds", pool=POOLS),
        metadata = config["outputs"]["output_dir"] + "expression_input/metadata.tagged.tsv.gz",
        wg1_psam = config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_psam"],
        exclude = lambda wildcards: config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/manual_selection/{cell_type}_PreQC_exclude_smf.txt".format(ancestry=wildcards.ancestry, cell_level=wildcards.cell_level, cell_type=wildcards.cell_type) if wildcards.qc == "PostQC" else []
    output:
        normalized_rds = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.Qced.Normalized.SCs.Rds",
        metadata = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.metadata.tsv.gz",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["quality_control"]["create_expression_matrices_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["quality_control"]["create_expression_matrices_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["quality_control"]["create_expression_matrices_time"]]
    threads: config["quality_control"]["create_expression_matrices_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/create_expression_matrices.R",
        input_dir = config["outputs"]["output_dir"] + "expression_input/pools/",
        exclude = lambda wildcards: "--exclude " + config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/manual_selection/{cell_type}_PreQC_exclude_smf.txt".format(ancestry=wildcards.ancestry, cell_level=wildcards.cell_level, cell_type=wildcards.cell_type) if wildcards.qc == "PostQC" else "",
        out_dir = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/",
    log: config["outputs"]["output_dir"] + "log/create_expression_matrices.{ancestry}.{cell_level}.{qc}.{cell_type}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --poolsheet {input.poolsheet} \
            --input_dir {params.input_dir} \
            --metadata {input.metadata} \
            --ancestry {wildcards.ancestry} \
            --cell_level {wildcards.cell_level} \
            --cell_type {wildcards.cell_type} \
            {params.exclude} \
            --out_dir {params.out_dir}
        """


# This functions takes a normalised expression Seurat object and applies to eQTL mapping filters and prepres the input files required for limix eQTL mapping.
# TODO: parameter aggregate_fun is not implemented
rule prepare_eqtl_input:
    input:
        seurat = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.Qced.Normalized.SCs.Rds",
        wg1_psam = config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_psam"],
        cell_annotation = config["inputs"]["cell_annotation_file"]
    output:
        smf = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.smf.txt",
        exp = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.Exp.txt",
        covariates = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.covariates.txt",
        pf = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.qtlInput.txt",
        cf = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.qtlInput.Pcs.txt",
        var = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.qtlInput.Pcs.var.txt",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["quality_control"]["prepare_eqtl_input_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["quality_control"]["prepare_eqtl_input_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["quality_control"]["prepare_eqtl_input_time"]]
    threads: config["quality_control"]["prepare_eqtl_input_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/prepare_eqtl_input.R",
        aggregate_fun = config["quality_control_extra"]["prepare_eqtl_input_aggregate_fun"],
        min_n_cells = config["quality_control_extra"]["prepare_eqtl_input_min_n_cells"],
        min_n_samples = config["quality_control_extra"]["prepare_eqtl_input_min_n_samples"],
        n_pcs = config["quality_control_extra"]["prepare_eqtl_input_n_pcs"],
        out_dir = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/",
    log: config["outputs"]["output_dir"] + "log/prepare_eqtl_input.{ancestry}.{cell_level}.{qc}.{cell_type}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --seurat {input.seurat} \
            --wg1_psam {input.wg1_psam} \
            --cell_annotation {input.cell_annotation} \
            --aggregate_fun {params.aggregate_fun} \
            --min_n_cells {params.min_n_cells} \
            --min_n_samples {params.min_n_samples} \
            --n_pcs {params.n_pcs} \
            --out_dir {params.out_dir} \
            --out_file {wildcards.cell_type}
        """


def get_kinship_file(wildcards):
    input_kinship = config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_kinship"].replace("{ancestry}", wildcards.ancestry)
    if os.path.exists(input_kinship):
        return input_kinship
    else:
        return config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}_sample.kinship"


# Generate a sample QC plot looking at cell count, cell fractions, expression PCs, and kinship PCs.
# Samples are filtered baed MAD thresholds selected by the user.
# TODO: consider checking the covariates file as well?
rule qtl_sample_qc:
    input:
        metadata = config["outputs"]["output_dir"] + "expression_input/metadata.tsv.gz",
        # covariates = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.covariates.txt",
        principal_components = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.qtlInput.Pcs.txt",
        kinship = get_kinship_file,
        smf = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.smf.txt",
    output:
        fig = report(config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/QC_figures/{cell_type}.{qc}.sample_qc.png", category="sample QC - {qc}", subcategory="{ancestry} - {cell_level} - {cell_type}", caption=config["inputs"]["repo_dir"] + "report_captions/sample_qc.rst"),
        sample_qc = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/manual_selection/{cell_type}_{qc}_sample_qc.txt",
        pass_qc_samples = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/manual_selection/{cell_type}_{qc}_exclude_smf.txt",
        thresholds = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/manual_selection/{cell_type}_{qc}_threshold_selection.txt"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["quality_control"]["qtl_sample_qc_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["quality_control"]["qtl_sample_qc_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["quality_control"]["qtl_sample_qc_time"]]
    threads: config["quality_control"]["qtl_sample_qc_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/qtl_sample_qc.py",
        donor_cell_threshold = "--donor_cell_threshold " + str(config["quality_control_extra"]["qtl_sample_qc_donor_cell_threshold"]) if config["quality_control_extra"]["qtl_sample_qc_donor_cell_threshold"] is not None else "",
        counts_threshold = "--counts_threshold " + str(config["quality_control_extra"]["qtl_sample_qc_counts_threshold"]) if config["quality_control_extra"]["qtl_sample_qc_counts_threshold"] is not None else "",
        kinship_threshold = "--kinship_threshold " + str(config["quality_control_extra"]["qtl_sample_qc_kinship_threshold"]) if config["quality_control_extra"]["qtl_sample_qc_kinship_threshold"] is not None else "",
        out = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/"
    log: config["outputs"]["output_dir"] + "log/qtl_sample_qc.{ancestry}.{cell_level}.{cell_type}.{qc}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --metadata {input.metadata} \
            --cell_level {wildcards.cell_level} \
            --cell_type {wildcards.cell_type} \
            --ancestry {wildcards.ancestry} \
            --qc {wildcards.qc} \
            {params.donor_cell_threshold} \
            --principal_components {input.principal_components} \
            {params.counts_threshold} \
            --kinship {input.kinship} \
            {params.kinship_threshold} \
            --smf {input.smf} \
            --out {params.out}
        """
