import os

CHROM = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']
configfile: "wp3_input_create.yaml"

image_folder = config["image_folder"]
top_dir = config["top_dir"]
wg1_folder = config["WG1_folder"]
wg2_folder = config["WG2_folder"]
wg3_folder = config["WG3_folder"]
singlet_assigned_seurat = config["wg1SingletsAssigned"]
wg1_psam = config["wg1_psam"]
wg1_wg2_combined_tag =  config["wg1_wg2_qc_taged"]
genotype_folder = wg3_folder+config["genotype_dir"]
genotype_merged_by_ancestry = config["ancestry_split_vcf"]
genotype_input_file = f"{wg1_folder}{genotype_merged_by_ancestry}EUR_imputed_hg38.vcf.gz"

scripts_folder = "/tools/WG3-pipeline-QTL/scripts/" ##Need to replace this.
unimputed_file = config["unimputed_folder"]
raw_genotype = f"{unimputed_file}"

expressionMetaD = wg3_folder+"input/AllMetaData.debug.txt"
sampleMappingFile = wg3_folder+'input/smf.txt'

kinshipFile = wg3_folder+"input/sample.kinship"
genotypeSampleInfo = genotype_folder+"EUR_imputed_hg38_stats.samples.gz"
cellAnnotationFile = config["cellAnnotation"]
inclusion_files = expand(genotype_folder+"EUR_imputed_hg38_inclusion_{chrom}.vars", chrom=CHROM)

annoFile = wg3_folder+'input/LimixAnnotationFile.txt'
chunkFile = wg3_folder+'input/ChunkingFile.txt'

if not os.path.exists(wg3_folder):
   os.makedirs(wg3_folder)

rule all:
    input:
        expand(genotype_folder+"EUR_imputed_hg38_varFiltered_chr{chrom}.bgen.metafile", chrom=CHROM), expressionMetaD, kinshipFile, chunkFile, annoFile, sampleMappingFile


rule Genotype_IO:
    input:
        genotype_input_file
    output:
        genotype_folder+"EUR_imputed_hg38_stats.vars.gz",
        genotype_folder+"EUR_imputed_hg38_stats.samples.gz",
        genotype_input_file+".tbi"
    shell:
        """
        singularity exec --bind {top_dir} {image_folder}wp3.simg tabix -p vcf  {genotype_input_file}
        singularity exec --bind {top_dir} {image_folder}wp3.simg java -Xmx5g -Xms5g -jar /tools/Genotype-IO-1.0.6-SNAPSHOT-jar-with-dependencies.jar -I VCF -i {genotype_input_file} -o {genotype_folder}EUR_imputed_hg38_stats
        gzip {genotype_folder}EUR_imputed_hg38_stats*
        """


rule ProcessOutputR:
    input:
        compressedVcf = genotype_folder+"EUR_imputed_hg38_stats.vars.gz",
        vcfIndex = genotype_input_file+".tbi"
    params:
        outD = genotype_folder
    output:
        temp(inclusion_files),
        outGenotypeLog = genotype_folder+"EUR_imputed_hg38_stats_filtered.vars.gz"
    shell:
        """
        singularity exec --bind {top_dir} {image_folder}wp3.simg Rscript {scripts_folder}process_output.R {input.compressedVcf} {params.outD}EUR_imputed_hg38_stats_filtered.vars {params.outD}EUR_imputed_hg38_inclusion_
        gzip {params.outD}EUR_imputed_hg38_stats_filtered.vars
        """


rule GenotypeHarmonizer:
    input:
        genotype_folder+"EUR_imputed_hg38_inclusion_{chrom}.vars"
    output:
        genotype_folder+"EUR_imputed_hg38_varFiltered_chr{chrom}.bgen"
    params:
        log=genotype_folder+"EUR_imputed_hg38_varFiltered_chr{chrom}.log",
        no_ext=genotype_folder+"EUR_imputed_hg38_varFiltered_chr{chrom}"
    shell:
        """
        singularity exec --bind {top_dir} {image_folder}wp3.simg java -Xmx10g -Xms10g -jar /tools/GenotypeHarmonizer-1.4.27-SNAPSHOT/GenotypeHarmonizer.jar -vf {input}  -i {genotype_input_file} -I VCF -O BGEN -o {params.no_ext} --genotypeField DS
        gzip {params.log}
        """


rule bgen_metadata_files:
    input:
        genotype_folder+"EUR_imputed_hg38_varFiltered_chr{chrom}.bgen"
    priority: 50
    output:
        temp(genotype_folder+"outputs/EUR_imputed_hg38_varFiltered_chr{chrom}.bgen.z"),
        temp(genotype_folder+"outputs/EUR_imputed_hg38_varFiltered_chr{chrom}.bgen.sample"),
        temp(genotype_folder+"outputs/EUR_imputed_hg38_varFiltered_chr{chrom}.bgen_master.txt"),
        genotype_folder+"EUR_imputed_hg38_varFiltered_chr{chrom}.bgen.metafile"
    shell:
        f"singularity exec --bind {top_dir} {image_folder}wp3.simg python {scripts_folder}make_z.py {input} {genotype_folder}"


rule kinship:
    input:
        raw_genotype+".pgen",
        raw_genotype+".psam",
        raw_genotype+".pvar"
    params:
        mainOut = wg3_folder+"genotypeR/raw_filtered"
    output:
         temp(wg3_folder+"genotypeR/raw_filtered.king")
    run:
        shell("mkdir -p {wg3_folder}genotypeR")
        shell("singularity exec --bind {top_dir} {image_folder}wp3.simg /tools/plink2 --pfile {raw_genotype} --maf 0.05 --hwe 1e-6 --make-pgen --out {wg3_folder}genotypeR/raw_filtered")
        shell("singularity exec --bind {top_dir} {image_folder}wp3.simg /tools/plink2 --pfile {wg3_folder}genotypeR/raw_filtered --indep-pairwise 250 50 0.2 --bad-ld --out {wg3_folder}genotypeR/raw_filtered")
        shell("singularity exec --bind {top_dir} {image_folder}wp3.simg /tools/plink2 --pfile {wg3_folder}genotypeR/raw_filtered --extract {wg3_folder}genotypeR/raw_filtered.prune.in --make-king square --out {params.mainOut}")


rule kinshipR:
    input:
        wg3_folder+"genotypeR/raw_filtered.king"
    priority: 50
    output:
        kinshipFile
    run:
        shell("singularity exec --bind {top_dir} {image_folder}wp3.simg Rscript {scripts_folder}kinship.R {input} {input}.id {output}")
        shell("rm -r {wg3_folder}genotypeR/")


rule derived_summarized_metadata:
    input:
        wg2_folder
    params:
        wg3_folder
    output:
        wg3_folder+wg1_wg2_combined_tag
    run:
        shell("singularity exec --bind {top_dir} {image_folder}/wp3.simg Rscript /tools/wg1-qc_filtering/QC_statistics_final.R "
            " --in_dir {input} "
            " --out_dir {params} ")


rule derive_expression_matrices:
    input:
        wg1 = wg1_folder+singlet_assigned_seurat,
        psam_wg1 = wg1_folder+wg1_psam,
        imputedSampleInfo = genotypeSampleInfo,
        cellInfo = cellAnnotationFile,
        summarizedMetaD = wg3_folder+wg1_wg2_combined_tag
    params:
        outD = wg3_folder+"input/L1/"
    output:
        temp(wg3_folder+"input/L1/tmpFiltered.Seurat.Rds"),
        wg3_folder+"input/AllMetaData.debug.txt",
        wg3_folder+"input/smf.txt"
    run:
        shell("mkdir -p {params.outD}")
        shell("singularity exec --bind {top_dir} {image_folder}/wp3.simg Rscript /tools/WG3-pipeline-QTL/scripts/Create_QCed_CT_pseudobulk_Azi_L1.R "
            " --wg1_data {input.wg1} "
            " --wg1_psam {input.psam_wg1} "
            " --metadata {input.summarizedMetaD} "
            " --wg1_imputed_genotypes {input.imputedSampleInfo} "
            " --cell_annotation {input.cellInfo} "
            " --out_dir {params.outD} ")


rule derived_feature_annotation_and_chunks:
    input:
        config["gtf_annotation_file"]
    params:
        nChunks = config["number_of_genes_in_chunk"],
        outFolder = wg3_folder+"input/"
    output:
        wg3_folder+'input/LimixAnnotationFile.txt',
        wg3_folder+'input/ChunkingFile.txt'
    run:
        shell("singularity exec --bind {top_dir} {image_folder}/wp3.simg Rscript {scripts_folder}createFeatureAnnotation.R "
            " --in_gtf {input} "
            " --n_genes {params.nChunks} "
            " --out_dir {params.outFolder} ")
