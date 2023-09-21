#!/usr/bin/env python
CHROMOSOMES = [str(chr) for chr in range(1, 23)]

# TODO: remove trailing / if from file paths.

rule all:
    input:
        expand(config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.metafile", chr=CHROMOSOMES),
        config["outputs"]["output_dir"] + "input/AllMetaData.debug.txt",
        config["outputs"]["output_dir"] + "input/sample.kinship",
        config["outputs"]["output_dir"] + "input/ChunkingFile.txt",
        config["outputs"]["output_dir"] + "input/LimixAnnotationFile.txt",
        config["outputs"]["output_dir"] + "input/smf.txt",
        expand(config["outputs"]["output_dir"] + "input/WindowsFilesChecks/chr_{chr}_windows_defined.txt", chr=CHROMOSOMES)


# TODO: make gzip explicit
rule Genotype_IO:
    input:
        vcf = config["inputs"]["wg1_folder"] + config["inputs_extra"]["relative_ancestry_split_path"] + "EUR_imputed_hg38.vcf.gz"
    output:
        compressed_vcf = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_stats.vars.gz",
        vcf_samples = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_stats.samples.gz"
    resources:
        java_mem = lambda wildcards, attempt: attempt * config["create_input"]["genotype_io_java_memory"],
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["create_input"]["genotype_io_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["create_input"]["genotype_io_memory"]
    threads: config["create_input"]["genotype_io_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        jar = "/tools/Genotype-IO-1.0.6-SNAPSHOT-jar-with-dependencies.jar",
        out = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_stats"
    log: config["outputs"]["output_dir"] + "logs/genotype_io.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} tabix -p vcf {input.vcf}
        singularity exec --bind {params.bind} {params.sif} java -Xmx{resources.java_mem}g -Xms{resources.java_mem}g -jar {params.jar} \
            -I VCF \
            -i {input.vcf} \
            -o {params.out}
        gzip {params.out}*
        """


rule ProcessOutputR:
    input:
        compressed_vcf = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_stats.vars.gz",
    output:
        vars = temp(expand(config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_inclusion_{chr}.vars", chr=CHROMOSOMES)),
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["create_input"]["process_output_r_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["create_input"]["process_output_r_memory"]
    threads: config["create_input"]["process_output_r_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/tools/WG3-pipeline-QTL/scripts/process_output.R",
        out = config["outputs"]["output_dir"] + "genotype_input"
    log: config["outputs"]["output_dir"] + "logs/process_output_r.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            {input.compressed_vcf} \
            {params.out}EUR_imputed_hg38_stats_filtered.vars \
            {params.out}EUR_imputed_hg38_inclusion_
        gzip {params.out}EUR_imputed_hg38_stats_filtered.vars
        """


# TODO: WTF waarom gzip de log file?
rule GenotypeHarmonizer:
    input:
        vcf = config["inputs"]["wg1_folder"] + config["inputs_extra"]["relative_ancestry_split_path"] + "EUR_imputed_hg38.vcf.gz",
        vars = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_inclusion_{chr}.vars"
    output:
        config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr{chr}.bgen"
    resources:
        java_mem = lambda wildcards, attempt: attempt * config["create_input"]["genotype_harmonizer_java_memory"],
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["create_input"]["genotype_harmonizer_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["create_input"]["genotype_harmonizer_memory"]
    threads: config["create_input"]["genotype_harmonizer_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        jar = "/tools/GenotypeHarmonizer-1.4.27-SNAPSHOT/GenotypeHarmonizer.jar",
        no_ext = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr{chr}"
    log: config["outputs"]["output_dir"] + "logs/genotype_harmonizer.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} java -Xmx{resources.java_mem}g -Xms{resources.java_mem}g -jar {params.jar} \
            -vf {input.vars} \
            -i {input.vcf} \
            -I VCF \
            -O BGEN \
            -o {params.no_ext} \
            --genotypeField DS
        gzip {log}
        """


rule bgen_metadata_files:
    priority: 50
    input:
        config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen"
    output:
        temp(config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.z"),
        temp(config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.sample"),
        temp(config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen_master.txt"),
        config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.metafile"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["create_input"]["bgen_metadata_files_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["create_input"]["bgen_metadata_files_memory"]
    threads: config["create_input"]["bgen_metadata_files_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/tools/WG3-pipeline-QTL/scripts/make_z.py",
        out = config["outputs"]["output_dir"] + "genotype_input"
    log: config["outputs"]["output_dir"] + "logs/bgen_metadata_files.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            {input} \
            {params.out}"
        """


# TODO check input?
rule kinship:
    input:
        pgen = config["inputs"]["unimputed_folder_folder"] + ".pgen",
        psam = config["inputs"]["unimputed_folder_folder"] + ".psam",
        pvar = config["inputs"]["unimputed_folder_folder"] + ".pvar"
    output:
         king = temp(config["outputs"]["output_dir"]+"genotypeR/raw_filtered.king")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["create_input"]["kinship_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["create_input"]["kinship_memory"]
    threads: config["create_input"]["kinship_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        maf = config["create_input_extra"]["kinship_maf"],
        hwe = config["create_input_extra"]["kinship_hwe"],
        out = config["outputs"]["output_dir"] + "genotypeR"
    log: config["outputs"]["output_dir"] + "logs/kinship.log"
    shell:
        """
        mkdir -p {params.out}
        singularity exec --bind {params.bind} {params.sif} /tools/plink2 \
            --threads {threads} \
            --pgen {input.pgen} \
            --pvar {input.pvar} \
            --psam {input.psam} \
            --maf {params.maf} \
            --hwe {params.hwe} \
            --make-pgen \
            --out {params.out}/raw_filtered

        singularity exec --bind {params.bind} {params.sif} /tools/plink2 \
            --pfile {params.out}/raw_filtered \
            --indep-pairwise 250 50 0.2 \
            --bad-ld \
            --out {params.out}/raw_filtered

        singularity exec --bind {params.bind} {params.sif} /tools/plink2 \
            --pfile {params.out}/raw_filtered \
            --extract {params.out}/raw_filtered.prune.in \
            --make-king square \
            --out {params.mainOut}/raw_filtered
        """


rule kinshipR:
    priority: 50
    input:
        king = config["outputs"]["output_dir"]+"genotypeR/raw_filtered.king",
        king_id = config["outputs"]["output_dir"]+"genotypeR/raw_filtered.king.id"
    output:
        kinship = config["outputs"]["output_dir"]+"input/sample.kinship"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["create_input"]["kinship_r_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["create_input"]["kinship_r_memory"]
    threads: config["create_input"]["kinship_r_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/tools/WG3-pipeline-QTL/scripts/kinship.R",
        out = config["outputs"]["output_dir"] + "genotypeR"
    log: config["outputs"]["output_dir"] + "logs/kinship_r.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            {input.king} \
            {input.king_id} \
            {output.kinship}
        
        rm -r {params.out}
        """


rule derived_summarized_metadata:
    input:
        in_dir = config["wg2_folder"]
    output:
        config["outputs"]["output_dir"] + "WG1_WG2_summary/qc_tag.rds"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["create_input"]["derived_summarized_metadata_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["create_input"]["derived_summarized_metadata_memory"]
    threads: config["create_input"]["derived_summarized_metadata_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/tools/wg1-qc_filtering/QC_statistics_final.R",
        out = config["outputs"]["output_dir"]
    log: config["outputs"]["output_dir"] + "logs/derived_summarized_metadata.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {script} \
            " --in_dir {input.in_dir} "
            " --out_dir {params.out}
        """


rule derive_expression_matrices:
    input:
        wg1 = config["inputs"]["wg1_folder"] + config["relative_wg1_singlets_assigned"],
        psam_wg1 = config["inputs"]["wg1_folder"] + config["inputs"]["relative_wg1_psam"],
        summarized_meta_dir = config["outputs"]["output_dir"] + "WG1_WG2_summary/qc_tag.rds",
        imputed_sample_info = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_stats.samples.gz",
        cell_info = config["inputs"]["cell_annotation"],
    output:
        temp(config["outputs"]["output_dir"] + "input/L1/tmpFiltered.Seurat.Rds"),
        config["outputs"]["output_dir"] + "input/AllMetaData.debug.txt",
        config["outputs"]["output_dir"] + "input/smf.txt"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["create_input"]["derive_expression_matrices_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["create_input"]["derive_expression_matrices_memory"]
    threads: config["create_input"]["derive_expression_matrices_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "input/L1/",
        script = "/tools/WG3-pipeline-QTL/scripts/Create_QCed_CT_pseudobulk_Azi_L1.R"
    log: config["outputs"]["output_dir"] + "logs/derive_expression_matrices.log"
    shell:
        """
        mkdir -p {params.out}
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --wg1_data {input.wg1} \
            --wg1_psam {input.psam_wg1} \
            --metadata {input.summarized_meta_dir} \
            --wg1_imputed_genotypes {input.imputed_sample_info} \
            --cell_annotation {input.cell_info} \
            --out_dir {params.out}
        """


rule derived_feature_annotation_and_chunks:
    input:
        gtf = config["ref"]["gtf_annotation_file"]
    output:
        config["outputs"]["output_dir"] + "input/LimixAnnotationFile.txt",
        config["outputs"]["output_dir"] + "input/ChunkingFile.txt"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["create_input"]["derived_feature_annotation_and_chunks_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["create_input"]["derived_feature_annotation_and_chunks_memory"]
    threads: config["create_input"]["derived_feature_annotation_and_chunks_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/tools/WG3-pipeline-QTL/scripts/createFeatureAnnotation.R",
        n_genes = config["create_input_extra"]["derived_feature_annotation_and_chunks_n_chuncks"],
        out = config["outputs"]["output_dir"] + "input/"
    log: config["outputs"]["output_dir"] + "logs/derived_feature_annotation_and_chunks.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --in_gtf {input.gtf} \
            --n_genes {params.n_genes} \
            --out_dir {params.out} \
            --autosomes_only
        """


rule defWindows:
    input:
        feature_file = config["outputs"]["output_dir"] + "input/LimixAnnotationFile.txt"
    output:
        config["outputs"]["output_dir"] + "input/WindowsFilesChecks/chr_{chr}_windows_defined.txt"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["create_input"]["def_windows_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["create_input"]["def_windows_memory"]
    threads: config["create_input"]["def_windows_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/tools/WG3-pipeline-QTL/scripts/define_windows.py",
        out = config["outputs"]["output_dir"] + "input/WindowsFiles",
        out_checks = config["outputs"]["output_dir"] + "input/WindowsFilesChecks"
    log: config["outputs"]["output_dir"] + "logs/def_windows.log"
    shell:
        """
        mkdir -p {params.out}
        mkdir -p {params.out_checks}
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            {input.feature_file} \
            {wildcards.chrom} \
            {params.out}
        touch {output}
        """
