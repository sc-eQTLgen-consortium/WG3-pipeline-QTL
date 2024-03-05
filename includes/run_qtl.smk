#!/usr/bin/env python
import os
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
            --out_dir {params.out}
        """


def parse_qtl_chunk(wildcards):
    match = re.match("([0-9]{1,2}|X|Y|MT)_([0-9]+)_([0-9]+)", wildcards.qtl_chunk)

    return match.group(1), match.group(2), match.group(3)


def get_qtl_bgen(wildcards):
    chr, _, _ = parse_qtl_chunk(wildcards)
    return config["outputs"]["output_dir"] + "genotype_input/{ancestry}/{ancestry}_imputed_hg38_varFiltered_chr{chr}.bgen".format(ancestry=wildcards.ancestry, chr=chr)


def get_qtl_kinship_file(wildcards):
    input_kinship = config["inputs"]["wg1_genotype_folder"] + config["inputs_extra"]["relative_wg1_kinship"].replace("{ancestry}", wildcards.ancestry)
    if os.path.exists(input_kinship):
        return input_kinship
    else:
        return config["outputs"]["output_dir"] + "kinship/{ancestry}/{ancestry}_sample.kinship"


def get_qtl_genomic_range(wilcards):
    chr, start, end = parse_qtl_chunk(wilcards)
    return "{chr}:{start}-{end}".format(chr=chr, start=start, end=end)


# Note: add -vf if you want to only test certain variants
rule run_qtl_mapping:
    priority: 10
    input:
        chunks = config["outputs"]["output_dir"] + "annot_input/{ancestry}/eQTLChunkingFile.txt",
        bgen = get_qtl_bgen,
        annotation_file = config["outputs"]["output_dir"] + "annot_input/{ancestry}/LimixAnnotationFile.txt",
        phenotype_file = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/PostQC/{cell_type}.qtlInput.txt",
        covariates_file = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/PostQC/{cell_type}.qtlInput.Pcs.txt",
        randomeff_files = get_qtl_kinship_file,
        sample_mapping_file = config["outputs"]["output_dir"] + "expression_input/{ancestry}/{cell_level}/{cell_type}/PostQC/{cell_type}.smf.txt"
    output:
        h5 = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/qtl_results_{qtl_chunk}.h5",
        # permutation = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/Permutation_information_{feature}.pickle.gz",
        # snp_qc_metrics = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/snp_qc_metrics_naContaining_feature_{feature}.txt.gz",
        snp_metadata = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/snp_metadata_{qtl_chunk}.txt.gz",
        feature_metadata = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/feature_metadata_{qtl_chunk}.txt.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["run_qtl_mapping_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["run_qtl_mapping_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["run_qtl"]["run_qtl_mapping_time"]]
    threads: config["run_qtl"]["run_qtl_mapping_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["limix_singularity_image"],
        script = "/limix_qtl/Limix_QTL/run_QTL_analysis_metaAnalysis.py",
        bgen = lambda wildcards, input: input.bgen.rstrip(".bgen"),
        output_directory = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/",
        window = config["run_qtl_extra"]["window_size"],
        genomic_range = get_qtl_genomic_range,
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
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --bgen {params.bgen} \
            --annotation_file {input.annotation_file} \
            --phenotype_file {input.phenotype_file} \
            --output_directory {params.output_directory} \
            --window {params.window} \
            --genomic_range {params.genomic_range} \
            --covariates_file {input.covariates_file} \
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
            --cis
        
        singularity exec --bind {params.bind} {params.sif} touch {output.h5}
        singularity exec --bind {params.bind} {params.sif} touch {output.snp_metadata}
        singularity exec --bind {params.bind} {params.sif} touch {output.feature_metadata}
        """


rule minimal_postprocess:
    input:
        h5 = expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/qtl_results_{qtl_chunk}.h5", qtl_chunk=QTL_CHUNKS, allow_missing=True),
        snp_metadata = expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/snp_metadata_{qtl_chunk}.txt.gz", qtl_chunk=QTL_CHUNKS, allow_missing=True),
        feature_metadata = expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/feature_metadata_{qtl_chunk}.txt.gz", qtl_chunk=QTL_CHUNKS, allow_missing=True)
    output:
        out = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl_results_all.txt.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["minimal_postprocess_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["minimal_postprocess_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["run_qtl"]["minimal_postprocess_time"]]
    threads: config["run_qtl"]["minimal_postprocess_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["limix_singularity_image"],
        script = "/limix_qtl/Limix_QTL/post_processing/minimal_postprocess.py",
        input_dir = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/",
        ouput_dir = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/",
        minimal_reporting_p = 1.0,
        minimal_reporting_featureP = 1.0
    log: config["outputs"]["output_dir"] + "log/minimal_postprocess.{ancestry}.{cell_level}.{cell_type}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --input_dir {params.input_dir} \
            --ouput_dir {params.ouput_dir} \
            --single_file_output \
            --write_compressed \
            --minimimal_reporting_p {params.minimal_reporting_p} \
            --minimal_reporting_featureP {params.minimal_reporting_featureP}
        """


rule top_feature:
    input:
        h5 = expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/qtl_results_{qtl_chunk}.h5", qtl_chunk=QTL_CHUNKS, allow_missing=True),
        snp_metadata = expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/snp_metadata_{qtl_chunk}.txt.gz", qtl_chunk=QTL_CHUNKS, allow_missing=True),
        feature_metadata = expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/feature_metadata_{qtl_chunk}.txt.gz", qtl_chunk=QTL_CHUNKS, allow_missing=True)
    output:
        out = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/top_qtl_results_all.txt.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["top_feature_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["run_qtl"]["top_feature_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["run_qtl"]["top_feature_time"]]
    threads: config["run_qtl"]["top_feature_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["limix_singularity_image"],
        script = "/limix_qtl/Limix_QTL/post_processing/minimal_postprocess.py",
        input_dir = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/",
        ouput_dir = config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/",
        minimal_reporting_p = 1.0,
        minimal_reporting_featureP = 1.0
    log: config["outputs"]["output_dir"] + "log/top_feature.{ancestry}.{cell_level}.{cell_type}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --input_dir {params.input_dir} \
            --ouput_dir {params.ouput_dir} \
            --single_file_output \
            --write_compressed \
            --minimimal_reporting_p {params.minimal_reporting_p} \
            --minimal_reporting_featureP {params.minimal_reporting_featureP} \
            --top_feature_based
        """


# TODO: remove 'qtl' folder after compression.
rule compress_qtl:
    input:
        h5 = expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/qtl_results_{qtl_chunk}.h5", zip, qtl_chunk=QTL_CHUNKS, allow_missing=True),
        # permutation = expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/Permutation_information_{feature}.pickle.gz", zip, qtl_chunk=QTL_CHUNKS, allow_missing=True),
        # snp_qc_metrics = expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/snp_qc_metrics_naContaining_feature_{feature}.txt.gz", zip, qtl_chunk=QTL_CHUNKS, allow_missing=True),
        snp_metadata = expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/snp_metadata_{qtl_chunk}.txt.gz", zip, qtl_chunk=QTL_CHUNKS, allow_missing=True),
        feature_metadata = expand(config["outputs"]["output_dir"] + "output/{ancestry}/{cell_level}/{cell_type}/qtl/feature_metadata_{qtl_chunk}.txt.gz",zip, qtl_chunk=QTL_CHUNKS, allow_missing=True),
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
