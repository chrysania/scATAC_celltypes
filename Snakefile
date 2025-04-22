# scATAC_celltypes snakefile

import pandas as pd

samples = pd.read_table("config.tsv").set_index("sample_name", drop=False).to_dict(orient='index')
sample_names = list(samples.keys())
tissue_set = sorted({v["sample_type"] for v in samples.values()})
tissue_sample_pairs = [(v["sample_type"], k) for k, v in samples.items()]

rule all:
    input:
#        expand("objects/{sample}_peaks.rds", sample=sample_names),
        expand("data/{tissue}/.download_complete", tissue=tissue_set)
        
rule download:
    input:
        "data/{tissue}/download.sh"
    output:
        touch("data/{tissue}/.download_complete")
    shell:
        """
        cd data/{wildcards.tissue}
        sh download.sh
        touch .download_complete
        """

rule count_barcodes:
    input:
        download_done = "data/{tissue}/.download_complete",
        fragments = "data/{tissue}/{sample}/fragments.tsv.gz"
    output:
        barcode_counts="data/{tissue}/{sample}/barcode_counts.tsv",
        barcodes="data/{tissue}/{sample}/barcodes_atac.txt"
    params:
        ncells=10000
    shell:
        """
        fragtk count \
            -f {input.fragments} \
            -o {output.barcode_counts} \
            -n {params.ncells} \
            > {output.barcodes}
        """

rule call_peaks:
    input:
        download_done = "data/{tissue}/.download_complete",
        fragments = "data/{tissue}/{sample}/fragments.tsv.gz"
    output:
        peaks="data/{tissue}/{sample}/peaks.bed"
    shell:
        """
        macs2 callpeak \
            -f BED --nomodel --shift -100 --extsize 200 --name {wildcards.sample} \
            -t {input.fragments} \
            --outdir data/{wildcards.tissue}/{wildcards.sample}/

        # cut narrowPeak to bed file
        cut -f1-3 data/{wildcards.tissue}/{wildcards.sample}/{wildcards.sample}_peaks.narrowPeak > {output}
        """

rule peak_matrix:
    input:
        download_done = "data/{tissue}/.download_complete",
        fragments=lambda wc: f"data/{samples[wc.sample]['sample_type']}/{wc.sample}/fragments.tsv.gz",
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
        done="data/{tissue}/.download_complete",
        frags=lambda wc: f"data/{samples[wc.sample]['sample_type']}/{wc.sample}/fragments.tsv.gz",
        barcodes=lambda wc: f"data/{samples[wc.sample]['sample_type']}/{wc.sample}/barcodes_atac.txt",
        peak_counts=lambda wc: f"data/{samples[wc.sample]['sample_type']}/{wc.sample}/peaks/",
        annotations="data/annotations.rds"
    output:
        object="objects/{sample}_peaks.rds"
    params:
        nCount_ATAC_above=lambda wc: samples[wc.sample]["nCount_ATAC_above"],
        nCount_ATAC_below=lambda wc: samples[wc.sample]["nCount_ATAC_below"],
        TSS_above=lambda wc: samples[wc.sample]["TSS_above"],
        nucleosome_signal=lambda wc: samples[wc.sample]["nucleosome_signal"]
    script:
        "code/build_object.R"

#rule combine_object:

# rule annotate_celltypes:
# todo: annotate celltypes for multiome datasets

# end : integrated_{tissue}.rds


