# scATAC_celltypes snakefile

import pandas as pd

samples = pd.read_table("config.tsv").set_index("sample_name", drop=False).to_dict(orient='index')
tissue_names = list(samples.keys())
tissue_sample_pairs = [(v["sample_type"], k) for k, v in samples.items()]

rule all:
    input:
        expand("objects/{sample}_peaks.rds", sample=[s for t, s in tissue_sample_pairs])

        
rule download:
    input: "data/{tissue}/download.sh"
    output: 
        "data/{tissue}/{sample}/fragments.tsv.gz",
        "data/{tissue}/{sample}/gex.mtx",
        "data/{tissue}/{sample}/genes.tsv",
        "data/{tissue}/{sample}/rna_cells.txt"
    shell:
        """
        cd data/{wildcards.tissue}
        sh download.sh
        """

rule count_barcodes:
    input: "data/{tissue}/{sample}/fragments.tsv.gz"
    output:
        barcode_counts="data/{tissue}/{sample}/barcode_counts.tsv",
        barcodes="data/{tissue}/{sample}/barcodes_atac.txt"
    params:
        ncells=10000
    shell:
        """
        fragtk count \
            -f {input} \
            -o {output.barcode_counts} \
            -n {params.ncells} \
            > {output.barcodes}
        """

rule call_peaks:
    input: 
        "data/{tissue}/{sample}/fragments.tsv.gz"
    output: 
        peaks="data/{tissue}/{sample}/peaks.bed"
    shell:
        """
        macs2 callpeak \
            -f BED --nomodel --shift -100 --extsize 200 --name {wildcards.sample} \
            -t {input} \
            --outdir data/{wildcards.tissue}/{wildcards.sample}/

        # cut narrowPeak to bed file
        cut -f1-3 data/{wildcards.tissue}/{wildcards.sample}/{wildcards.sample}_peaks.narrowPeak > {output}
        """

rule peak_matrix:
    input:
        frags="data/{tissue}/{sample}/fragments.tsv.gz",
        regions="data/{tissue}/{sample}/peaks.bed",
        barcodes="data/{tissue}/{sample}/barcodes_atac.txt"
    output:
        directory("data/{tissue}/{sample}/peaks/")
    shell:
        """
        fragtk matrix \
            -f {input.frags} \
            -c {input.barcodes} \
            -o {output} \
            -b {input.regions} \
            --pic
        """

rule build_object:
    input:
        frags=lambda wildcards: f"data/{samples[wildcards.sample]['sample_type']}/{wildcards.sample}/fragments.tsv.gz",
        barcodes=lambda wildcards: f"data/{samples[wildcards.sample]['sample_type']}/{wildcards.sample}/barcodes_atac.txt",
        peak_counts=lambda wildcards: f"data/{samples[wildcards.sample]['sample_type']}/{wildcards.sample}/peaks/",
        annotations="data/annotations.rds"
    output:
        object="objects/{sample}_peaks.rds"
    params:
        nCount_ATAC_above=lambda wildcards: samples[wildcards.sample]["nCount_ATAC_above"],
        nCount_ATAC_below=lambda wildcards: samples[wildcards.sample]["nCount_ATAC_below"],
        TSS_above=lambda wildcards: samples[wildcards.sample]["TSS_above"],
        nucleosome_signal=lambda wildcards: samples[wildcards.sample]["nucleosome_signal"]
    script:
        "code/build_object.R"

#rule combine_object:

# rule annotate_celltypes:
# todo: annotate celltypes for multiome datasets

# end : integrated_{tissue}.rds


