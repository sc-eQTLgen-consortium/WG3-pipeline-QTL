#!/usr/bin/env python
import re
import os

# TODO: remove trailing / if from file paths.

CELLTYPES = ["L1/AST", "L1/END", "L1/EX", "L1/IN", "L1/MIC", "L1/OLI", "L1/OPC", "L1/PER"]

chunk_chr, chunk_start, chunk_end=[], [], []
with open(config["inputs"]["wg3_folder"] + "input/ChunkingFile.txt") as fp:
    for line in fp:
        re_match=re.match(r"([A-Za-z0-9]+):(\d+)-(\d+)", line.strip())
        chunk_chr.append(re_match[1])
        chunk_start.append(re_match[2])
        chunk_end.append(re_match[3])

rule all:
    input:
        expand(config["outputs"]["output_dir"] + "{ct}/top_qtl_results_all.txt.gz", ct=CELLTYPES),
        expand(config["outputs"]["output_dir"] + "{ct}/qtl.h5.tgz", ct=CELLTYPES),
        expand(config["outputs"]["output_dir"] + "{ct}/qtl.annotation.tgz", ct=CELLTYPES),
        expand(config["outputs"]["output_dir"] + "{ct}/qtl.permutations.tgz", ct=CELLTYPES),
        config["outputs"]["output_dir"] + "LDMatrices.tgz"
    output:
        touch(expand(config["outputs"]["output_dir"] + "{ct}/done.txt", ct=CELLTYPES))


# TODO: -rf = kf, was that correct?
rule run_qtl_mapping:
    priority: 10
    input:
        af = config["inputs"]["wg3_folder"] + "input/LimixAnnotationFile.txt",
        cf = config["inputs"]["wg3_folder"] + "input/{ct}.qtlInput.Pcs.txt",
        pf = config["inputs"]["wg3_folder"] + "input/{ct}.qtlInput.txt",
        rf = config["inputs"]["wg3_folder"] + "input/sample.kinship",
        smf = config["inputs"]["wg3_folder"] + "input/smf.txt"
    output:
        touch(config["outputs"]["output_dir"] + "{ct}/qtl/{chr}_{start}_{end}.finished")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["run_qtl_mapping_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["run_qtl_mapping_memory"]
    threads: config["run_qtl"]["run_qtl_mapping_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/limix_qtl/Limix_QTL/run_QTL_analysis_metaAnalysis.py",
        gen = config["inputs"]["wg3_folder"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}",
        od = config["outputs"]["output_dir"] + "{ct}/qtl/",
        np = config["run_qtl_extra"]["n_permutations"],
        maf = config["run_qtl_extra"]["maf"],
        w = config["run_qtl_extra"]["window_size"],
        hwe = config["run_qtl_extra"]["hwe_cutoff"],
        rs = config["run_qtl_extra"]["relatedness_score"],
    log: config["outputs"]["output_dir"] + "logs/run_qtl_mapping.{ct}.chr_{chr}_{start}_{end}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --bgen {params.gen} \
            -af {input.af} \
            -cf {input.cf} \
            -pf {input.pf} \
            -rf {input.rf} \
            -smf {input.smf} \
            -od {params.od} \
            -gr {wildcards.chr}:{wildcards.start}-{wildcards.end} \
            -np {params.np}\
            -maf {params.maf}\
            -c \
            -gm gaussnorm \
            -w {params.w} \
            -hwe {params.hwe} \
            -rs {params.rs}
        """


# TODO: does this zip work? looks weird.
rule top_feature:
    input:
        expand(config["outputs"]["output_dir"]+ "{ct}/qtl/{chr}_{start}_{end}.finished", zip, chr=chunk_chr, start=chunk_start, end=chunk_end, allow_missing=True)
    output:
        temp(config["outputs"]["output_dir"] + "{ct}/top_qtl_results_all.txt")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["top_feature_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["top_feature_memory"]
    threads: config["run_qtl"]["top_feature_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/limix_qtl/Limix_QTL/post_processing/minimal_postprocess.py",
        in_dir = config["outputs"]["output_dir"] + "{ct}/qtl/",
        out = config["outputs"]["output_dir"] + "{ct}/"
    log: config["outputs"]["output_dir"] + "logs/top_feature.{ct}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            -id {params.in_dir} \
            -od {params.out} \
            -tfb \
            -sfo
        """


rule compress_qtl:
    input:
        topQtl = config["outputs"]["output_dir"] + "{ct}/top_qtl_results_all.txt",
        qtlChunks = expand(config["outputs"]["output_dir"]+ "{ct}/qtl/{chr}_{start}_{end}.finished", zip, chr=chunk_chr, start=chunk_start, end=chunk_end, allow_missing=True)
    output:
        h5 = config["outputs"]["output_dir"] + "{ct}/qtl.h5.tgz",
        anno = config["outputs"]["output_dir"] + "{ct}/qtl.annotation.tgz",
        perm = config["outputs"]["output_dir"] + "{ct}/qtl.permutations.tgz",
        gz_top_qtl = config["outputs"]["output_dir"] + "{ct}/top_qtl_results_all.txt.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["compress_qtl_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["compress_qtl_memory"]
    threads: config["run_qtl"]["compress_qtl_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        dir = config["outputs"]["output_dir"] + "{ct}/qtl"
    log: config["outputs"]["output_dir"] + "logs/compress_qtl.{ct}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} tar -czf {output.h5} {params.dir}/*.h5 && find {params.dir} | grep '\.pickle.gz$' > {params.dir}/files.txt &&  tar cf {output.perm} -T {params.dir}/files.txt && tar -cf {output.anno} {params.dir}/*.txt.gz && gzip {input.topQtl}
        """


rule make_temporary_files:
    input:
        bgen_file = config["inputs"]["wg3_folder"]+ "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen"
    output:
        temp(config["inputs"]["wg3_folder"]+ "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.z"),
        temp(config["inputs"]["wg3_folder"]+ "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.sample"),
        temp(config["inputs"]["wg3_folder"]+ "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen_master.txt")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["make_temporary_files_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["make_temporary_files_memory"]
    threads: config["run_qtl"]["make_temporary_files_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/tools/WG3-pipeline-QTL/scripts/make_z.py",
        in_dir = config["inputs"]["wg3_folder"] + "genotype_input"
    log: config["outputs"]["output_dir"] + "logs/make_temporary_files.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            {input.bgen_file} \
            {params.in_dir}
        """


rule create_bdose_file_by_chr:
    input:
        config["inputs"]["wg3_folder"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen",
        config["inputs"]["wg3_folder"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.bgi",
        config["inputs"]["wg3_folder"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.z",
        config["inputs"]["wg3_folder"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.sample",
        in_files = config["inputs"]["wg3_folder"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen_master.txt"
    output:
        temp(config["inputs"]["wg3_folder"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.bdose"),
        temp(config["inputs"]["wg3_folder"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.bdose.bdose.tmp0"),
        temp(config["inputs"]["wg3_folder"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.bdose.meta.tmp0")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["create_bdose_file_by_chr_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["create_bdose_file_by_chr_memory"]
    threads: config["run_qtl"]["create_bdose_file_by_chr_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        ldstore = "/tools/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64"
    log: config["outputs"]["output_dir"] + "logs/create_bdose_file_by_chr.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} {params.ldstore} \
            --in-files {input.in_files} \
            --write-bdose \
            --bdose-version 1.1
        """

rule LD_per_window:
    input:
        bdose_file = config["inputs"]["wg3_folder"] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.bdose",
        window = config["inputs"]["wg3_folder"] + "input/WindowsFiles/chr_{chr}_window_{num}_{start}_{end}.pkl"
    output:
        config["inputs"]["wg3_folder"] + "input/LDMatrices/chr_{chr}_window_{num}_{start}_{end}_means.pkl.gz",
        config["inputs"]["wg3_folder"] + "input/LDMatrices/chr_{chr}_window_{num}_{start}_{end}_low_dim.pkl.gz",
        config["inputs"]["wg3_folder"] + "input/LDMatrices/chr_{chr}_window_{num}_{start}_{end}_components.pkl.gz",
        touch(config["inputs"]["wg3_folder"] + "input/LDMatrices/chr_{chr}_window_{num}_{start}_{end}.pkl.gz")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["ld_per_window_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["ld_per_window_memory"]
    threads: config["run_qtl"]["ld_per_window_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/tools/WG3-pipeline-QTL/scripts/process_LD_windows.py",
        out_dir = config["inputs"]["wg3_folder"] + "input",
        out = config["inputs"]["wg3_folder"] + "input/LDMatrices/chr_{chr}_window_{num}_{start}_{end}.pkl"
    log: config["outputs"]["output_dir"] + "logs/ld_per_window.chr_{chr}_window_{num}_{start}_{end}.log"
    shell:
        """
        mkdir -p {params.out_dir}/LDMatrices
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            {params.out_dir} \
            {input.bdose_file} \
            {input.window} \
            {wildcards.chr} \
            {params.out}
        """


rule compress_ld:
    input:
        ldWindows = expand(config["inputs"]["wg3_folder"] + "input/LDMatrices/{winF}.gz", winF=os.listdir(config["inputs"]["wg3_folder"] + "input/WindowsFiles"))
    output:
        config["outputs"]["output_dir"] + "LDMatrices.tgz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["compress_ld_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["compress_ld_memory"]
    threads: config["run_qtl"]["compress_ld_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["inputs"]["wg3_folder"] + "input/LDMatrices"
    log: config["outputs"]["output_dir"] + "logs/compress_ld.log"
    shell: 
        """
        tar -czf {output} {params.out}"
        """
