#!/usr/bin/env python
import os

######################################
############# ANNOTATION #############
######################################

rule feature_annotation:
    input:
        gtf = config["refs"]["gtf_annotation_file"]
    output:
        feature_file = config["outputs"]["output_dir"] + "annot_input/{ancestry}/LimixAnnotationFile.txt"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["feature_annotation_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["feature_annotation_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["prepare_input"]["feature_annotation_time"]]
    threads: config["prepare_input"]["feature_annotation_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/create_feature_file.R",
        feature_name = config["prepare_input_extra"]["feature_annotation_feature_name"],
        biotype_flag = config["prepare_input_extra"]["feature_annotation_biotype_flag"],
        autosomes_only = "--autosomes_only " if config["prepare_input_extra"]["feature_annotation_autosomes_only"] else "",
        out = config["outputs"]["output_dir"] + "annot_input/{ancestry}/"
    log: config["outputs"]["output_dir"] + "log/feature_annotation.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --in_gtf {input.gtf} \
            --feature_name {params.feature_name} \
            --biotype_flag {params.biotype_flag} \
            {params.autosomes_only} \
            --out_dir {params.out} > {log} 2>&1
        """

#####################################
############# GENOTYPES #############
#####################################


rule get_genotype_samples:
    input:
        wg1_metadata = config["inputs"]["wg1_demultiplex_folder"] + config["inputs_extra"]["relative_wg1_metadata"]
    output:
        samples = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_individuals.txt"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["get_genotype_samples_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["get_genotype_samples_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["prepare_input"]["get_genotype_samples_time"]]
    threads: config["prepare_input"]["get_genotype_samples_threads"]
    params:
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"],
        script = config["inputs"]["repo_dir"] + "scripts/get_genotype_samples.py",
        individual_aggregate = config["settings_extra"]["individual_aggregate"],
        output_dir = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/",
    log: config["outputs"]["output_dir"] + "log/get_genotype_samples.{ancestry}.log"
    shell:
        """
         singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --wg1_metadata {input.wg1_metadata} \
            --ancestry {wildcards.ancestry} \
            --individual_aggregate {params.individual_aggregate} \
            --output_dir {params.output_dir} > {log} 2>&1
        """


rule filter_vcf:
    input:
        vcf = config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_imputed_genotype_vcf"],
        samples = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_individuals.txt"
    output:
        vcf = temp(config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38.vcf.gz"),
        index1 = temp(config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38.vcf.gz.csi"),
        index2 = temp(config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38.vcf.gz.tbi")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["filter_vcf_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["filter_vcf_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["prepare_input"]["filter_vcf_time"]]
    threads: config["prepare_input"]["filter_vcf_threads"]
    params:
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"],
        script = config["inputs"]["repo_dir"] + "scripts/extract_samples_from_metadata.py",
        ouput_dir = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/",
    log: config["outputs"]["output_dir"] + "log/filter_vcf.{ancestry}.log"
    shell:
        """
         singularity exec --bind {params.bind} {params.sif} bcftools view -S {input.samples} {input.vcf} --force-samples -Oz -o {output.vcf}
         singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
         singularity exec --bind {params.bind} {params.sif} tabix -p vcf {output.vcf}
        """


def get_input_vcf(wildcards):
    if config["settings"]["filter_samples"]:
        return config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38.vcf.gz"
    else:
        return config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_imputed_genotype_vcf"]

def get_input_vcf_index(wildcards):
    vcf = get_input_vcf(wildcards=wildcards)
    return vcf + ".tbi"


rule genotype_io:
    input:
        vcf = get_input_vcf,
        index = get_input_vcf_index
    output:
        vars = temp(config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_stats.vars.gz"),
        samples = temp(config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_stats.samples.gz")
    resources:
        java_mem_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["genotype_io_memory"] * config["prepare_input"]["genotype_io_threads"] - config["settings_extra"]["java_memory_buffer"],
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["genotype_io_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["genotype_io_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["prepare_input"]["genotype_io_time"]]
    threads: config["prepare_input"]["genotype_io_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        jar = "/opt/Genotype-IO-1.0.6-SNAPSHOT/GenotypeIO.jar ", # Used to be: /tools/Genotype-IO-1.0.6-SNAPSHOT-jar-with-dependencies.jar
        out = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_stats"
    log: config["outputs"]["output_dir"] + "log/genotype_io.{ancestry}.log"
    shell:
        """
         singularity exec --bind {params.bind} {params.sif} java -Xmx{resources.java_mem_gb}g -Xms{resources.java_mem_gb}g -jar {params.jar} \
            --inputType VCF \
            --input {input.vcf} \
            --output {params.out}
         singularity exec --bind {params.bind} {params.sif} gzip {params.out}.vars
         singularity exec --bind {params.bind} {params.sif} gzip {params.out}.samples
        """


rule filter_variants:
    input:
        vars = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_stats.vars.gz"
    output:
        vars = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_stats_filtered.vars.gz",
        vars_per_chr = temp(expand(config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_inclusion_{chr}.vars", chr=CHROMOSOMES, allow_missing=True))
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["filter_variants_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["filter_variants_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["prepare_input"]["filter_variants_time"]]
    threads: config["prepare_input"]["filter_variants_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/filter_variants.py",
        maf = config["prepare_input_extra"]["filter_variants_maf"],
        r2 = config["prepare_input_extra"]["filter_variants_r2"],
        call = config["prepare_input_extra"]["filter_variants_call"],
        out = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/"
    log: config["outputs"]["output_dir"] + "log/filter_variants.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --input {input.vars} \
            --maf {params.maf} \
            --r2 {params.r2} \
            --call {params.call} \
            --out {params.out} > {log} 2>&1
        """


# Added a check here to see that the variant is not empty. GenotypeHarmonizer doesn't care.
rule genotype_harmonizer:
    input:
        vcf = get_input_vcf,
        index = get_input_vcf_index,
        vars = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_inclusion_{chr}.vars"
    output:
        bgen = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.bgen",
        bgi = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.bgen.bgi",
        log = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.log"
    resources:
        java_mem_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["genotype_harmonizer_memory"] * config["prepare_input"]["genotype_harmonizer_threads"] - config["settings_extra"]["java_memory_buffer"],
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["genotype_harmonizer_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["genotype_harmonizer_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["prepare_input"]["genotype_harmonizer_time"]]
    threads: config["prepare_input"]["genotype_harmonizer_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        jar = "/opt/GenotypeHarmonizer-1.4.27-SNAPSHOT/GenotypeHarmonizer.jar", # Used to be: /tools/GenotypeHarmonizer-1.4.27-SNAPSHOT/GenotypeHarmonizer.jar
        output = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}"
    log: config["outputs"]["output_dir"] + "log/genotype_harmonizer.{ancestry}.chr_{chr}.log"
    shell:
        """
        if [[ "$( singularity exec --bind {params.bind} {params.sif} cat {input.vars} | wc -l)" -eq "0" ]]; 
        then
           echo "Error, the variants file is empty"
           exit
        fi
        
         singularity exec --bind {params.bind} {params.sif} java -Xmx{resources.java_mem_gb}g -Xms{resources.java_mem_gb}g -jar {params.jar} \
            --variantFilterList {input.vars} \
            --input {input.vcf} \
            --inputType VCF \
            --outputType BGEN \
            --output {params.output} \
            --genotypeField DS
        """


rule bgen_metadata_files:
    priority: 50
    input:
        bgen = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.bgen",
        index = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.bgen.bgi"
    output:
        z = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.z",
        sample = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.sample",
        master = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}_master.txt",
        bgen_metafile = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.bgen.metafile"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["bgen_metadata_files_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["bgen_metadata_files_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["prepare_input"]["bgen_metadata_files_time"]]
    threads: config["prepare_input"]["bgen_metadata_files_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        python = "/opt/conda/envs/py38/bin/python3.8",  # Used to be: python
        script = config["inputs"]["repo_dir"] + "scripts/make_z.py",
        out = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}"
    log: config["outputs"]["output_dir"] + "log/bgen_metadata_files.{ancestry}.chr_{chr}.log"
    shell:
        """
         singularity exec --bind {params.bind} {params.sif} {params.python} {params.script} \
            --in_filepath {input.bgen} \
            --out {params.out} > {log} 2>&1
        """


rule kinship:
    input:
        vcf = config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_het_filtered"],
        samples = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_samples.txt" if config["settings"]["filter_samples"] else []
    output:
        subset_pgen = temp(config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}_subset.pgen"),
        subset_pvar = temp(config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}_subset.pvar"),
        subset_psam = temp(config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}_subset.psam"),
        subset_log = config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}_subset.log",
        prune_in = temp(config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}_subset_pruning.prune.in"),
        prune_out = temp(config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}_subset_pruning.prune.out"),
        prune_log = config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}_subset_pruning.log",
        pruned_log = config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}_subset_pruned.log",
        king = temp(config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}_subset_pruned.king"),
        king_id = temp(config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}_subset_pruned.king.id"),
        kinship = config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}_sample.kinship"
    resources:
        plink_mem_mb = lambda wildcards, attempt: (attempt * config["prepare_input"]["kinship_memory"] * config["prepare_input"]["kinship_threads"] - config["settings_extra"]["plink_memory_buffer"]) * 1000,
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["kinship_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["kinship_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["prepare_input"]["kinship_time"]]
    threads: config["prepare_input"]["kinship_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        plink = "plink2", # Used to be: /tools/plink2
        genome_build = config["settings"]["genome_build"],
        maf = config["prepare_input_extra"]["kinship_maf"],
        hwe = config["prepare_input_extra"]["kinship_hwe"],
        out_subset = config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}_subset",
        indep_window_size = config["prepare_input_extra"]["kinship_indep_window_size"],
        indep_step_size = config["prepare_input_extra"]["kinship_indep_step_size"],
        pairwise_r2_threshold = config["prepare_input_extra"]["kinship_pairwise_r2_threshold"],
        out_pruning = config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}_subset_pruning",
        keep = lambda wildcards, input: "--keep " + input.samples if config["settings"]["filter_samples"] else "",
        out_pruned = config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}_subset_pruned",
        script = config["inputs"]["repo_dir"] + "scripts/kinship.R"
    log: config["outputs"]["output_dir"] + "log/kinship.{ancestry}.log"
    shell:
        """        
         singularity exec --bind {params.bind} {params.sif} {params.plink} \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --vcf {input.vcf} \
            --split-par {params.genome_build} \
            --maf {params.maf} \
            --hwe {params.hwe} \
            {params.keep} \
            --make-pgen \
            --out {params.out_subset}

         singularity exec --bind {params.bind} {params.sif} {params.plink} \
            --memory {resources.plink_mem_mb} \
            --pgen {output.subset_pgen} \
            --pvar {output.subset_pvar} \
            --psam {output.subset_psam} \
            --indep-pairwise {params.indep_window_size} {params.indep_step_size} {params.pairwise_r2_threshold} \
            --bad-ld \
            --out {params.out_pruning}

         singularity exec --bind {params.bind} {params.sif} {params.plink} \
            --memory {resources.plink_mem_mb} \
            --pgen {output.subset_pgen} \
            --pvar {output.subset_pvar} \
            --psam {output.subset_psam} \
            --extract {output.prune_in} \
            --make-king square \
            --out {params.out_pruned}
            
         singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --king {output.king} \
            --king_id {output.king_id} \
            --out {output.kinship}
        """


def get_kinship_file(wildcards):
    input_kinship = config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_kinship"].replace("{ancestry}", wildcards.ancestry)
    if os.path.exists(input_kinship):
        return input_kinship
    else:
        return config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}_sample.kinship"


rule kinship_pca:
    input:
        kinship = get_kinship_file,
        wg1_psam = config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_psam"],
    output:
        pca = config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}.kinship.Pcs.txt",
        pca_var = config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}.kinship.Pcs.var.txt",
        pca_plot = report(config["outputs"]["output_dir"] + "kinship/{ancestry}/QC_figures/{ancestry}.kinship.Pcs.png", category="kinship QC", subcategory="{ancestry}", caption=config["inputs"]["repo_dir"] + "report_captions/kinship_pca.rst"),
        pca_var_plot = report(config["outputs"]["output_dir"] + "kinship/{ancestry}/QC_figures/{ancestry}.kinship.Pcs.var.png", category="kinship QC", subcategory="{ancestry}", caption=config["inputs"]["repo_dir"] + "report_captions/kinship_pca_var.rst"),
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["kinship_pca_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["kinship_pca_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["prepare_input"]["kinship_pca_time"]]
    threads: config["prepare_input"]["kinship_pca_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/kinship_pca.R",
        data_out = config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}.kinship.",
        plot_out = config["outputs"]["output_dir"] + "kinship/{ancestry}/QC_figures/{ancestry}.kinship."
    log: config["outputs"]["output_dir"] + "log/kinship_pca.{ancestry}.log"
    shell:
        """        
         singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --kinship {input.kinship} \
            --wg1_psam {input.wg1_psam} \
            --data_out {params.data_out} \
            --plot_out {params.plot_out} > {log} 2>&1
        """


######################################
############# EXPRESSION #############
######################################


# WG1: imputed VCF file + metadata file (pool, barcode, singlet / doublet, assignment) per pool
# WG2: metadata file (pool, barcode, cell annotation)
# L2 -> L1 mapping is added.
#
# Calculate nfeatures_RNA, ncounts_RNA, complexity (nfeatures_RNA) / ncounts_RNA, percent.mt, percent.rb per pool
rule create_seurat:
    input:
        poolsheet = config["inputs"]["poolsheet_path"],
        cell_annotation = config["inputs"]["cell_annotation_file"],
        wg1_metadata = config["inputs"]["wg1_demultiplex_folder"] + config["inputs_extra"]["relative_wg1_metadata"],
        wg1_psam = config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_psam"],
        wg2_metadata = config["inputs"]["wg2_folder"] + config["inputs_extra"]["relative_wg2_metadata"],
        rb_genes = config["refs"]["ref_dir"] + config["refs_extra"]["ribosomal_genes"],
        mt_genes = config["refs"]["ref_dir"] + config["refs_extra"]["mitochondrial_genes"]
    output:
        seurat = config["outputs"]["output_dir"] + "expression_input/pools/{pool}.rds",
        full_metadata = config["outputs"]["output_dir"] + "expression_input/pools/{pool}.full.metadata.tsv.gz",
        qc_metrics = config["outputs"]["output_dir"] + "expression_input/pools/{pool}.qc_metrics.tsv.gz",
        filter_stats = config["outputs"]["output_dir"] + "expression_input/pools/{pool}.filter.stats.tsv",
        metadata = config["outputs"]["output_dir"] + "expression_input/pools/{pool}.metadata.tsv.gz",
        tmpdir = temp(directory(config["outputs"]["output_dir"] + "tmp/create_seurat/{pool}/")),
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["create_seurat_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["create_seurat_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["prepare_input"]["create_seurat_time"]],
    threads: config["prepare_input"]["create_seurat_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/create_seurat.R",
        wg2_pairing = "--wg2_pairing " + config["refs"]["ref_dir"] + config["refs_extra"]["wg2_pairing"] if config["refs_extra"]["wg2_pairing"] is not None and os.path.exists(config["refs"]["ref_dir"] + config["refs_extra"]["wg2_pairing"]) else "",
        out = config["outputs"]["output_dir"] + "expression_input/pools/"
    log: config["outputs"]["output_dir"] + "log/create_seurat.{pool}.log"
    shell:
        """
         export TMPDIR={output.tmpdir}
         mkdir -p $TMPDIR
        
         singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --poolsheet {input.poolsheet} \
            --pool {wildcards.pool} \
            --cell_annotation {input.cell_annotation} \
            --wg1_metadata {input.wg1_metadata} \
            --wg2_metadata {input.wg2_metadata} \
            {params.wg2_pairing} \
            --wg1_psam {input.wg1_psam} \
            --rb_genes {input.rb_genes} \
            --mt_genes {input.mt_genes} \
            --out_dir {params.out} > {log} 2>&1
        
        singularity exec --bind {params.bind} {params.sif} touch {output.metadata} {output.seurat}
            
         unset TMPDIR
        """


# In this rule, we want to read this meta data file and calculate MAD stuff. This has to be over all pools together!
# A tag is added to the metadata to mark if the barcode passed QC or not.
# TODO: parameter md_vars and downsampling is not implemented
# TODO: do per cell type ?
rule barcode_qc:
    input:
        poolsheet = config["inputs"]["poolsheet_path"],
        metadata = expand(config["outputs"]["output_dir"] + "expression_input/pools/{pool}.metadata.tsv.gz", pool=POOLS)
    output:
        metadata = config["outputs"]["output_dir"] + "expression_input/metadata.tsv.gz",
        tagged = config["outputs"]["output_dir"] + "expression_input/metadata.tagged.tsv.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["barcode_qc_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["barcode_qc_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["prepare_input"]["barcode_qc_time"]]
    threads: config["prepare_input"]["barcode_qc_threads"]
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
            --out_dir {params.out_dir} > {log} 2>&1
        """

rule barcode_qc_plots:
    input:
        metadata = config["outputs"]["output_dir"] + "expression_input/metadata.tsv.gz"
    output:
        tmpdir = temp(directory(config["outputs"]["output_dir"] + "tmp/barcode_qc_plots/")),
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
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["barcode_qc_plots_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["barcode_qc_plots_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["prepare_input"]["barcode_qc_plots_time"]]
    threads: config["prepare_input"]["barcode_qc_plots_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/barcode_qc_plots.R",
        out = config["outputs"]["output_dir"] + "expression_input/QC_figures/"
    log: config["outputs"]["output_dir"] + "log/barcode_qc_plots.log"
    shell:
        """
        export TMPDIR={output.tmpdir}
        mkdir -p $TMPDIR

        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --metadata {input.metadata} \
            --out {params.out} > {log} 2>&1

        singularity exec --bind {params.bind} {params.sif} touch {output.done}

        unset TMPDIR
        """


# For each cell type load the pool seurat object, filter on singlet, barcode QC tag, cell type, ancestry (optional), treatment and combine into one seurat object.
# Then apply PF - log - PF: calculate total number of UMI over all samples (meanSampleSum).
# Optional: sample are removed baed on rule 'sample_qc'.
rule create_expression_matrices:
    input:
        poolsheet = config["inputs"]["poolsheet_path"],
        seurat = expand(config["outputs"]["output_dir"] + "expression_input/pools/{pool}.rds", pool=POOLS),
        filter_stats = expand(config["outputs"]["output_dir"] + "expression_input/pools/{pool}.filter.stats.tsv", pool=POOLS),
        metadata = config["outputs"]["output_dir"] + "expression_input/metadata.tagged.tsv.gz",
        wg1_psam = config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_psam"],
        exclude = lambda wildcards: config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/manual_selection/{cell_type}_PreQC_exclude_smf.txt".format(ancestry=wildcards.ancestry, cell_level=wildcards.cell_level, cell_type=wildcards.cell_type) if wildcards.qc == "PostQC" else []
    output:
        filter_stats = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.filter.stats.tsv.gz",
        normalized_rds = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.Qced.Normalized.SCs.Rds",
        metadata = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/{qc}/{cell_type}.metadata.tsv.gz",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["create_expression_matrices_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["create_expression_matrices_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["prepare_input"]["create_expression_matrices_time"]]
    threads: config["prepare_input"]["create_expression_matrices_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/create_expression_matrices.R",
        input_dir = config["outputs"]["output_dir"] + "expression_input/pools/",
        exclude = lambda wildcards, input: "--exclude " + input.exclude if input.exclude else "",
        individual_aggregate = config["settings_extra"]["individual_aggregate"],
        sample_aggregate = config["settings_extra"]["sample_aggregate"],
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
            --individual_aggregate {params.individual_aggregate} \
            --sample_aggregate {params.sample_aggregate} \
            --out_dir {params.out_dir} > {log} 2>&1
        """
