#!/usr/bin/env python
import os

rule run_qtl_mapping:
    priority: 10
    input:
        af = config["outputs"]["output_dir" ] + "input/LimixAnnotationFile.txt",
        cf = config["outputs"]["output_dir" ] + "input/{ct}.qtlInput.Pcs.txt",
        pf = config["outputs"]["output_dir" ] + "input/{ct}.qtlInput.txt",
        kf = config["outputs"]["output_dir" ] + "input/sample.kinship",
        smf = config["outputs"]["output_dir" ] + "input/l1/smf.txt"
    output:
        touch(config["outputs"]["output_dir"] + "output/{ct}/qtl/chr_{chr}_{start}_{end}.finished")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["run_qtl_mapping_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["run_qtl_mapping_memory"]
    threads: config["run_qtl"]["run_qtl_mapping_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["limix_singularity_image"],
        script = "/limix_qtl/Limix_QTL/run_QTL_analysis_metaAnalysis.py",
        gen = config["outputs"]["output_dir" ] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}",
        od = config["outputs"]["output_dir"] + "output/{ct}/qtl/",
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
            -rf {input.kf} \
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


rule top_feature:
    input:
        expand(config["outputs"]["output_dir"]+ "output/{ct}/qtl/chr_{chr}_{start}_{end}.finished", ct=CELLTYPES, chr=CHUNK_CHR, start=CHUNK_START, end=CHUNK_END)
    output:
        temp(config["outputs"]["output_dir"] + "output/{ct}/top_qtl_results_all.txt")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["top_feature_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["top_feature_memory"]
    threads: config["run_qtl"]["top_feature_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["limix_singularity_image"],
        script = "/limix_qtl/Limix_QTL/post_processing/minimal_postprocess.py",
        in_dir = config["outputs"]["output_dir"] + "output/{ct}/qtl/",
        out = config["outputs"]["output_dir"] + "output/{ct}/"
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
        topQtl = config["outputs"]["output_dir"] + "output/{ct}/top_qtl_results_all.txt",
        qtlChunks = expand(config["outputs"]["output_dir"]+ "output/{ct}/qtl/chr_{chr}_{start}_{end}.finished", ct=CELLTYPES, chr=CHUNK_CHR, start=CHUNK_START, end=CHUNK_END)
    output:
        h5 = config["outputs"]["output_dir"] + "output/{ct}/qtl.h5.tgz",
        anno = config["outputs"]["output_dir"] + "output/{ct}/qtl.annotation.tgz",
        perm = config["outputs"]["output_dir"] + "output/{ct}/qtl.permutations.tgz",
        gz_top_qtl = config["outputs"]["output_dir"] + "output/{ct}/top_qtl_results_all.txt.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["compress_qtl_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["compress_qtl_memory"]
    threads: config["run_qtl"]["compress_qtl_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        dir = config["outputs"]["output_dir"] + "output/{ct}/qtl"
    log: config["outputs"]["output_dir"] + "logs/compress_qtl.{ct}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} tar -czf {output.h5} {params.dir}/*.h5 && find {params.dir} | grep '\.pickle.gz$' > {params.dir}/files.txt &&  tar cf {output.perm} -T {params.dir}/files.txt && tar -cf {output.anno} {params.dir}/*.txt.gz && gzip {input.topQtl}
        """


rule create_bdose_file_by_chr:
    input:
        config["outputs"]["output_dir" ] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen",
        config["outputs"]["output_dir" ] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.bgi",
        config["outputs"]["output_dir" ] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.z",
        config["outputs"]["output_dir" ] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.sample",
        in_files = config["outputs"]["output_dir" ] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen_master.txt"
    output:
        temp(config["outputs"]["output_dir" ] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.bdose"),
        temp(config["outputs"]["output_dir" ] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.bdose.bdose.tmp0"),
        temp(config["outputs"]["output_dir" ] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.bdose.meta.tmp0")
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
        bdose_file = config["outputs"]["output_dir" ] + "genotype_input/EUR_imputed_hg38_varFiltered_chr_{chr}.bgen.bdose",
        window = config["outputs"]["output_dir" ] + "input/WindowsFiles/chr_{chr}_window_{num}_{start}_{end}.pkl"
    output:
        low_dim = config["outputs"]["output_dir"] + "input/LDMatrices/chr_{chr}_window_{num}_{start}_{end}_low_dim.pkl.gz",
        components = config["outputs"]["output_dir" ] + "input/LDMatrices/chr_{chr}_window_{num}_{start}_{end}_components.pkl.gz",
        means = config["outputs"]["output_dir"] + "input/LDMatrices/chr_{chr}_window_{num}_{start}_{end}_means.pkl.gz",
        data = touch(config["outputs"]["output_dir" ] + "input/LDMatrices/chr_{chr}_window_{num}_{start}_{end}.pkl.gz")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["ld_per_window_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["ld_per_window_memory"]
    threads: config["run_qtl"]["ld_per_window_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/tools/WG3-pipeline-QTL/scripts/process_LD_windows.py",
        out_dir = config["outputs"]["output_dir" ] + "input",
        out = config["outputs"]["output_dir" ] + "input/LDMatrices/chr_{chr}_window_{num}_{start}_{end}.pkl"
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


# TODO: the os.listdir is needed since we don't know the numbering of the windows. However this seems a bit risky since it may include
#  junk files in the compressed file that are not intended to be there. The commented fix gives a 'No values given for wildcard 'num'' error
#  since I don't know the window numbering when defining the required input_files for rule all.
rule compress_ld:
    input:
        dWindows = expand(config["outputs"]["output_dir"] + "input/LDMatrices/{winF}.gz", winF=os.listdir(config["outputs"]["output_dir"] + "input/WindowsFiles"))
        # low_dim = expand(config["outputs"]["output_dir"] + "input/LDMatrices/chr_{chr}_window_{num}_{start}_{end}_low_dim.pkl.gz", ct=CELLTYPES, chr=CHUNK_CHR, start=CHUNK_START, end=CHUNK_END),
        # components = expand(config["outputs"]["output_dir" ] + "input/LDMatrices/chr_{chr}_window_{num}_{start}_{end}_components.pkl.gz", ct=CELLTYPES, chr=CHUNK_CHR, start=CHUNK_START, end=CHUNK_END),
        # means = expand(config["outputs"]["output_dir"] + "input/LDMatrices/chr_{chr}_window_{num}_{start}_{end}_means.pkl.gz", ct=CELLTYPES, chr=CHUNK_CHR, start=CHUNK_START, end=CHUNK_END),
        # data = expand(config["outputs"]["output_dir" ] + "input/LDMatrices/chr_{chr}_window_{num}_{start}_{end}.pkl.gz", ct=CELLTYPES, chr=CHUNK_CHR, start=CHUNK_START, end=CHUNK_END),
    output:
        config["outputs"]["output_dir"] + "output/LDMatrices.tgz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["compress_ld_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["compress_ld_memory"]
    threads: config["run_qtl"]["compress_ld_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir" ] + "LDMatrices.tgz"
    log: config["outputs"]["output_dir"] + "logs/compress_ld.log"
    shell:
        """
        tar -czf {params.out} {input.low_dim} {input.components} {input.means} {input.data}
        """
