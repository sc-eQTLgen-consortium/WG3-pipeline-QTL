#!/usr/bin/env python
import os

rule genotype_io:
    input:
        vcf = config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_imputed_genotype_vcf"],
        index = config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_imputed_genotype_vcf"] + ".tbi",
    output:
        vars = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_stats.vars.gz",
        samples = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_stats.samples.gz"
    resources:
        java_mem_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["genotype_io_memory"] * config["prepare_input"]["genotype_io_threads"] - config["settings_extra"]["java_memory_buffer"],
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["genotype_io_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["genotype_io_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["prepare_input"]["genotype_io_time"]]
    threads: config["prepare_input"]["genotype_io_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        jar = "/tools/Genotype-IO-1.0.6-SNAPSHOT-jar-with-dependencies.jar",
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


rule process_output_r:
    input:
        vars = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_stats.vars.gz"
    output:
        vars = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_stats_filtered.vars.gz",
        vars_per_chr = temp(expand(config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_inclusion_{chr}.vars", chr=CHROMOSOMES, allow_missing=True))
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["process_output_r_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["process_output_r_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["prepare_input"]["process_output_r_time"]]
    threads: config["prepare_input"]["process_output_r_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/process_output.R",
        maf = config["prepare_input_extra"]["process_output_r_maf"],
        r2 = config["prepare_input_extra"]["process_output_r_r2"],
        call = config["prepare_input_extra"]["process_output_r_call"],
        out = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/"
    log: config["outputs"]["output_dir"] + "log/process_output_r.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --input {input.vars} \
            --maf {params.maf} \
            --r2 {params.r2} \
            --call {params.call} \
            --out {params.out}
        """


rule genotype_harmonizer:
    input:
        vcf = config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_imputed_genotype_vcf"],
#        index = config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_imputed_genotype_vcf"] + ".csi",
        vars = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_inclusion_{chr}.vars"
    output:
        log = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.log",
        bgen = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.bgen",
        bgi = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.bgen.bgi"
    resources:
        java_mem_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["genotype_harmonizer_memory"] * config["prepare_input"]["genotype_harmonizer_threads"] - config["settings_extra"]["java_memory_buffer"],
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["genotype_harmonizer_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["genotype_harmonizer_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["prepare_input"]["genotype_harmonizer_time"]]
    threads: config["prepare_input"]["genotype_harmonizer_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        jar = "/tools/GenotypeHarmonizer-1.4.27-SNAPSHOT/GenotypeHarmonizer.jar",
        output = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}"
    log: config["outputs"]["output_dir"] + "log/genotype_harmonizer.{ancestry}.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} java -Xmx{resources.java_mem_gb}g -Xms{resources.java_mem_gb}g -jar {params.jar} \
            --variantFilterList {input.vars} \
            --input {input.vcf} \
            --inputType VCF \
            --outputType BGEN \
            --output {params.output} \
            --genotypeField DS
        """


rule kinship:
    input:
        vcf = config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_het_filtered"]
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
        genome_build = config["settings"]["genome_build"],
        maf = config["prepare_input_extra"]["kinship_maf"],
        hwe = config["prepare_input_extra"]["kinship_hwe"],
        out_subset = config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}_subset",
        indep_window_size = config["prepare_input_extra"]["kinship_indep_window_size"],
        indep_step_size = config["prepare_input_extra"]["kinship_indep_step_size"],
        pairwise_r2_threshold = config["prepare_input_extra"]["kinship_pairwise_r2_threshold"],
        out_pruning = config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}_subset_pruning",
        out_pruned = config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}_subset_pruned",
        script = config["inputs"]["repo_dir"] + "scripts/kinship.R",
    log: config["outputs"]["output_dir"] + "log/kinship.{ancestry}.log"
    shell:
        """        
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --vcf {input.vcf} \
            --split-par {params.genome_build} \
            --maf {params.maf} \
            --hwe {params.hwe} \
            --make-pgen \
            --out {params.out_subset}

        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --pgen {output.subset_pgen} \
            --pvar {output.subset_pvar} \
            --psam {output.subset_psam} \
            --indep-pairwise {params.indep_window_size} {params.indep_step_size} {params.pairwise_r2_threshold} \
            --bad-ld \
            --out {params.out_pruning}

        singularity exec --bind {params.bind} {params.sif} plink2 \
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
        metadata = config["outputs"]["output_dir"] + "expression_input/pools/{pool}.metadata.tsv.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["create_seurat_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["create_seurat_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["prepare_input"]["create_seurat_time"]]
    threads: config["prepare_input"]["create_seurat_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/create_seurat.R",
        wg2_pairing = "--wg2_pairing " + config["refs"]["ref_dir"] + config["refs_extra"]["wg2_pairing"] if os.path.exists(config["refs"]["ref_dir"] + config["refs_extra"]["wg2_pairing"]) else "",
        out = config["outputs"]["output_dir"] + "expression_input/pools/"
    log: config["outputs"]["output_dir"] + "log/create_seurat.{pool}.log"
    shell:
        """
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
            --out_dir {params.out}
        """


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
            --out_dir {params.out}
        """
