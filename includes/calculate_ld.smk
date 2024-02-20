#!/usr/bin/env python
import re

LD_CHUNK_PATTERN = "chr_([0-9]{1,2}|X|Y|MT)_window_([0-9]+)_([0-9]+)_([0-9]+)"


rule bgen_metadata_files:
    priority: 50
    input:
        bgen = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.bgen",
        index = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.bgen.bgi"
    output:
        z = temp(config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.z"),
        sample = temp(config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.sample"),
        master = temp(config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}_master.txt")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["calculate_ld"]["bgen_metadata_files_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["calculate_ld"]["bgen_metadata_files_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["calculate_ld"]["bgen_metadata_files_time"]]
    threads: config["calculate_ld"]["bgen_metadata_files_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/make_z.py",
        out = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}"
    log: config["outputs"]["output_dir"] + "log/bgen_metadata_files.{ancestry}.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --in_filepath {input.bgen} \
            --out {params.out}
        """


rule create_bdose_file_by_chr:
    input:
        bgen = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.bgen",
        index = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.bgen.bgi",
        z = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.z",
        sample = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.sample",
        master = config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}_master.txt"
    output:
        bdose = temp(config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.bdose"),
        bdose_bdose = temp(config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.bdose.bdose.tmp0"),
        bdose_meta = temp(config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.bdose.meta.tmp0")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["calculate_ld"]["create_bdose_file_by_chr_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["calculate_ld"]["create_bdose_file_by_chr_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["calculate_ld"]["create_bdose_file_by_chr_time"]]
    threads: config["calculate_ld"]["create_bdose_file_by_chr_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        ldstore = "/tools/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64"
    log: config["outputs"]["output_dir"] + "log/create_bdose_file_by_chr.{ancestry}.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} {params.ldstore} \
            --in-files {input.master} \
            --write-bdose \
            --bdose-version 1.1
        """


rule ld_chunks:
    input:
        feature_file = config["outputs"]["output_dir"] + "annot_input/{ancestry}/LimixAnnotationFile.txt"
    output:
        chunks = config["outputs"]["output_dir"] + "annot_input/{ancestry}/LDChunkingFile.txt"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["calculate_ld"]["ld_chunks_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["calculate_ld"]["ld_chunks_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["calculate_ld"]["ld_chunks_time"]]
    threads: config["calculate_ld"]["ld_chunks_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/derived_ld_chunks.py",
        out = config["outputs"]["output_dir"] + "annot_input/{ancestry}/",
        window_size = config["calculate_ld_extra"]["ld_chunks_window_size"],
        gene_window = config["calculate_ld_extra"]["ld_chunks_gene_window"],
    log: config["outputs"]["output_dir"] + "log/ld_chunks.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --feature_file {input.feature_file} \
            --window_size {params.window_size} \
            --gene_window {params.gene_window} \
            --out {params.out}
        """


def parse_ld_chunk(wildcards):
    match = re.match(LD_CHUNK_PATTERN, wildcards.ld_chunk)
    return match.group(1), match.group(2), match.group(3), match.group(4)


def get_ld_bdose(wildcards):
    chr, _, _, _ = parse_ld_chunk(wildcards)
    return config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.bdose".format(ancestry=wildcards.ancestry, chr=chr)


rule LD_per_window:
    input:
        chunks = config["outputs"]["output_dir"] + "annot_input/{ancestry}/LDChunkingFile.txt",
        bdose = get_ld_bdose
    output:
        low_dim = config["outputs"]["output_dir"] + "output/{ancestry}/LDMatrices/chr{ld_chr}/{ld_chunk}_low_dim.pkl.gz",
        components = config["outputs"]["output_dir"] + "output/{ancestry}/LDMatrices/chr{ld_chr}/{ld_chunk}_components.pkl.gz",
        means = config["outputs"]["output_dir"] + "output/{ancestry}/LDMatrices/chr{ld_chr}/{ld_chunk}_means.pkl.gz",
        data = config["outputs"]["output_dir"] + "output/{ancestry}/LDMatrices/chr{ld_chr}/{ld_chunk}.pkl.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["calculate_ld"]["ld_per_window_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["calculate_ld"]["ld_per_window_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["calculate_ld"]["ld_per_window_time"]]
    threads: config["calculate_ld"]["ld_per_window_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/process_LD_windows.py",
        out = config["outputs"]["output_dir" ] + "output/{ancestry}/LDMatrices/chr{ld_chr}/"
    log: config["outputs"]["output_dir"] + "log/ld_per_window.{ancestry}.chr{ld_chr}.{ld_chunk}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --bdose {input.bdose} \
            --outfile {wildcards.ld_chunk} \
            --out {params.out}
        """



rule compress_ld:
    input:
        low_dim = expand(config["outputs"]["output_dir"] + "output/{ancestry}/LDMatrices/chr{ld_chr}/{ld_chunk}_low_dim.pkl.gz", zip, ld_chr=[re.match(LD_CHUNK_PATTERN, ld_chunk).group(1) for ld_chunk in LD_CHUNKS], ld_chunk=LD_CHUNKS, allow_missing=True),
        components = expand(config["outputs"]["output_dir"] + "output/{ancestry}/LDMatrices/chr{ld_chr}/{ld_chunk}_components.pkl.gz", zip, ld_chr=[re.match(LD_CHUNK_PATTERN, ld_chunk).group(1) for ld_chunk in LD_CHUNKS], ld_chunk=LD_CHUNKS, allow_missing=True),
        means = expand(config["outputs"]["output_dir"] + "output/{ancestry}/LDMatrices/chr{ld_chr}/{ld_chunk}_means.pkl.gz", zip, ld_chr=[re.match(LD_CHUNK_PATTERN, ld_chunk).group(1) for ld_chunk in LD_CHUNKS], ld_chunk=LD_CHUNKS, allow_missing=True),
        data = expand(config["outputs"]["output_dir"] + "output/{ancestry}/LDMatrices/chr{ld_chr}/{ld_chunk}.pkl.gz", zip, ld_chr=[re.match(LD_CHUNK_PATTERN, ld_chunk).group(1) for ld_chunk in LD_CHUNKS], ld_chunk=LD_CHUNKS, allow_missing=True),
    output:
        out = config["outputs"]["output_dir"] + "output/{ancestry}/LDMatrices/LDMatrices.tgz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["calculate_ld"]["compress_ld_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["calculate_ld"]["compress_ld_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["calculate_ld"]["compress_ld_time"]]
    threads: config["calculate_ld"]["compress_ld_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        merge_list = config["outputs"]["output_dir"] + "output/{ancestry}/LDMatrices/merge_list.txt"
    log: config["outputs"]["output_dir"] + "log/compress_ld.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} echo {input.low_dim} | sed 's/ /\\n/g' > {params.merge_list}
        singularity exec --bind {params.bind} {params.sif} echo {input.components} | sed 's/ /\\n/g' >> {params.merge_list}
        singularity exec --bind {params.bind} {params.sif} echo {input.means} | sed 's/ /\\n/g' >> {params.merge_list}
        singularity exec --bind {params.bind} {params.sif} echo {input.data} | sed 's/ /\\n/g' >> {params.merge_list}
        singularity exec --bind {params.bind} {params.sif} tar -czf {output.out} -T {params.merge_list}
        """