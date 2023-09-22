#!/usr/bin/env python


rule Genotype_IO:
    input:
        vcf = config["inputs"]["wg1_folder"] + config["inputs_extra"]["relative_ancestry_split_path"] + "EUR_imputed_hg38.vcf.gz",
        index = config["inputs"]["wg1_folder"] + config["inputs_extra"]["relative_ancestry_split_path"] + "EUR_imputed_hg38.vcf.gz.tbi"
    output:
        vars_gz = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_stats.vars.gz",
        samples_gz = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_stats.samples.gz"
    resources:
        java_mem = lambda wildcards, attempt: attempt * config["prepare_input"]["genotype_io_java_memory"],
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["genotype_io_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["genotype_io_memory"]
    threads: config["prepare_input"]["genotype_io_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        jar = "/tools/Genotype-IO-1.0.6-SNAPSHOT-jar-with-dependencies.jar",
        out = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_stats"
    log: config["outputs"]["output_dir"] + "logs/genotype_io.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} java -Xmx{resources.java_mem}g -Xms{resources.java_mem}g -jar {params.jar} \
            -I VCF \
            -i {input.vcf} \
            -o {params.out}
        gzip {params.out}.vars
        gzip {params.out}.samples
        """


rule ProcessOutputR:
    input:
        vars = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_stats.vars.gz",
    output:
        vars = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_stats_filtered.vars",
        vars_per_chr = temp(expand(config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_inclusion_chr_{chr}.vars", chr=CHROMOSOMES)),
        vars_gz = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_stats_filtered.vars.gz",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["process_output_r_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["process_output_r_memory"]
    threads: config["prepare_input"]["process_output_r_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/tools/WG3-pipeline-QTL/scripts/process_output.R",
        out = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_inclusion_chr_"
    log: config["outputs"]["output_dir"] + "logs/process_output_r.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            {input.vars} \
            {output.vars} \
            {params.out}
        gzip {output.vars}
        """


rule GenotypeHarmonizer:
    input:
        vcf = config["inputs"]["wg1_folder"] + config["inputs_extra"]["relative_ancestry_split_path"] + "EUR_imputed_hg38.vcf.gz",
        index = config["inputs"]["wg1_folder"] + config["inputs_extra"]["relative_ancestry_split_path"] + "EUR_imputed_hg38.vcf.gz.tbi",
        vars = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_inclusion_chr_{chr}.vars"
    output:
        bgen = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen",
        bgi = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.bgi"
    resources:
        java_mem = lambda wildcards, attempt: attempt * config["prepare_input"]["genotype_harmonizer_java_memory"],
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["genotype_harmonizer_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["genotype_harmonizer_memory"]
    threads: config["prepare_input"]["genotype_harmonizer_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        jar = "/tools/GenotypeHarmonizer-1.4.27-SNAPSHOT/GenotypeHarmonizer.jar",
        no_ext = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}"
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
        """


rule bgen_metadata_files:
    priority: 50
    input:
        bgen = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen",
        bgi = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.bgi"
    output:
        z = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.z",
        sample = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.sample",
        master = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen_master.txt",
        metafile = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.metafile"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["bgen_metadata_files_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["bgen_metadata_files_memory"]
    threads: config["prepare_input"]["bgen_metadata_files_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/tools/WG3-pipeline-QTL/scripts/make_z.py",
        out = config["outputs"]["output_dir"] + "genotype_input"
    log: config["outputs"]["output_dir"] + "logs/bgen_metadata_files.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            {input.bgen} \
            {params.out}"
        """


rule kinship:
    input:
        pgen = config["inputs"]["unimputed_folder"] + ".pgen",
        psam = config["inputs"]["unimputed_folder"] + ".psam",
        pvar = config["inputs"]["unimputed_folder"] + ".pvar"
    output:
        king = temp(config["outputs"]["output_dir" ] + "genotypeR/raw_filtered.king"),
        king_id = temp(config["outputs"]["output_dir"] + "genotypeR/raw_filtered.king.id")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["kinship_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["kinship_memory"]
    threads: config["prepare_input"]["kinship_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        maf = config["create_input_extra"]["kinship_maf"],
        hwe = config["create_input_extra"]["kinship_hwe"],
        out = config["outputs"]["output_dir"] + "genotypeR/raw_filtered"
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
            --out {params.out}

        singularity exec --bind {params.bind} {params.sif} /tools/plink2 \
            --pfile {params.out} \
            --indep-pairwise 250 50 0.2 \
            --bad-ld \
            --out {params.out}

        singularity exec --bind {params.bind} {params.sif} /tools/plink2 \
            --pfile {params.out} \
            --extract {params.out}.prune.in \
            --make-king square \
            --out {params.out}
        """


rule kinshipR:
    priority: 50
    input:
        king = config["outputs"]["output_dir" ] + "genotypeR/raw_filtered.king",
        king_id = config["outputs"]["output_dir" ] + "genotypeR/raw_filtered.king.id"
    output:
        kinship = config["outputs"]["output_dir" ] + "input/sample.kinship"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["kinship_r_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["kinship_r_memory"]
    threads: config["prepare_input"]["kinship_r_threads"]
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
        """


# TODO: switch back the rds files as required input. Can't do this right now since my dataset does not have these files.
rule derived_summarized_metadata:
    input:
        # metadata_reduced_data = config["inputs"]["wg2_folder"] + "step4_reduce/metadata.reduced_data.RDS",
        # reduced_data = config["inputs"]["wg2_folder"] + "step4_reduce/reduced_data.RDS",
        in_dir = config["inputs"]["wg2_folder"]
    output:
        qc_tag = config["outputs"]["output_dir"] + "WG1_WG2_summary/qc_tag.rds"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["derived_summarized_metadata_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["derived_summarized_metadata_memory"]
    threads: config["prepare_input"]["derived_summarized_metadata_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/tools/wg1-qc_filtering/QC_statistics_final.R",
        in_dir = config["inputs"]["wg2_folder"],
        out = config["outputs"]["output_dir"]
    log: config["outputs"]["output_dir"] + "logs/derived_summarized_metadata.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {script} \
            " --in_dir {params.in_dir} "
            " --out_dir {params.out}
        """


rule derive_expression_matrices:
    input:
        wg1_data = config["inputs"]["wg1_folder"] + config["inputs_extra"]["relative_wg1_singlets_assigned"],
        wg1_psam = config["inputs"]["wg1_folder"] + config["inputs_extra"]["relative_wg1_psam"],
        metadata = config["outputs"]["output_dir"] + "WG1_WG2_summary/qc_tag.rds",
        wg1_imputed_genotypes = config["outputs"]["output_dir"] + "genotype_input/EUR_imputed_hg38_stats.samples.gz",
        cell_annotation = config["inputs"]["cell_annotation"],
    output:
        tmp_filtered_rds = temp(config["outputs"]["output_dir"] + "input/L1/tmpFiltered.Seurat.Rds"),
        normalized_rds = expand(config["outputs"]["output_dir"] + "input/{ct}.Qced.Normalized.SCs.Rds",ct=CELLTYPES),
        exp = expand(config["outputs"]["output_dir"] + "input/{ct}.Exp.txt", ct=CELLTYPES),
        covariates = expand(config["outputs"]["output_dir"] + "input/{ct}.covariates.txt",ct=CELLTYPES),
        pf = expand(config["outputs"]["output_dir"] + "input/{ct}.qtlInput.txt", ct=CELLTYPES),
        cf = expand(config["outputs"]["output_dir"] + "input/{ct}.qtlInput.Pcs.txt",ct=CELLTYPES),
        metadata = config["outputs"]["output_dir"] + "input/L1/AllMetaData.debug.txt",
        smf = config["outputs"]["output_dir"] + "input/L1/smf.txt"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["derive_expression_matrices_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["derive_expression_matrices_memory"]
    threads: config["prepare_input"]["derive_expression_matrices_threads"]
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
            --wg1_data {input.wg1_data} \
            --wg1_psam {input.wg1_psam} \
            --metadata {input.metadata} \
            --wg1_imputed_genotypes {input.wg1_imputed_genotypes} \
            --cell_annotation {input.cell_annotation} \
            --out_dir {params.out}
        """


rule derived_feature_annotation_and_chunks:
    input:
        gtf = config["refs"]["gtf_annotation_file"]
    output:
        annot_file = config["outputs"]["output_dir"] + "input/LimixAnnotationFile.txt",
        chunk_file = config["outputs"]["output_dir"] + "input/ChunkingFile.txt"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["derived_feature_annotation_and_chunks_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["derived_feature_annotation_and_chunks_memory"]
    threads: config["prepare_input"]["derived_feature_annotation_and_chunks_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/tools/WG3-pipeline-QTL/scripts/createFeatureAnnotation.R",
        n_genes = config["create_input_extra"]["derived_feature_annotation_and_chunks_n_genes"],
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


# WindowsFilesChecks are not the real files outputted but just tmp files to trick
# snakemake
rule defWindows:
    input:
        feature_file = config["outputs"]["output_dir"] + "input/LimixAnnotationFile.txt"
    output:
        windows = config["outputs"]["output_dir"] + "input/WindowsFilesChecks/chr_{chr}_windows_defined.txt",
        # windows = expand(config["outputs"]["output_dir"] + "input/WindowsFiles/chr_{chr}_window_{num}_{start}_{end}.pkl", ct=CELLTYPES, chr=CHUNK_CHR, start=CHUNK_START, end=CHUNK_END)
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["def_windows_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["prepare_input"]["def_windows_memory"]
    threads: config["prepare_input"]["def_windows_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/tools/WG3-pipeline-QTL/scripts/define_windows.py",
        out = config["outputs"]["output_dir"] + "input/WindowsFiles",
        out_checks = config["outputs"]["output_dir"] + "input/WindowsFilesChecks"
    log: config["outputs"]["output_dir"] + "logs/def_windows.chr_{chr}.log"
    shell:
        """
        mkdir -p {params.out}
        mkdir -p {params.out_checks}
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            {input.feature_file} \
            {wildcards.chrom} \
            {params.out}
        touch {output.windows}
        """
