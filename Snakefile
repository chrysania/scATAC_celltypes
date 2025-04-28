# scATAC_celltypes snakefile

import pandas as pd

samples = pd.read_table("config.tsv").set_index("sample_name", drop=False).to_dict(orient='index')
sample_names = list(samples.keys())
tissue_set = sorted({v["sample_type"] for v in samples.values()})
tissue_sample_pairs = [(v["sample_type"], k) for k, v in samples.items()]

from statistics import mean

tissue_clustering_resolution = {
    tissue: mean([
        v["clustering_resolution"] for k, v in samples.items()
        if v["sample_type"] == tissue
    ])
    for tissue in tissue_set
}

rule all:
    input:
#        expand("data/{tissue}/.download_complete", tissue=tissue_set)
        expand("objects/{sample}_peaks.rds", sample=sample_names),
        expand("data/{tissue}/{sample}/peaks/", zip, tissue=[x[0] for x in tissue_sample_pairs], sample=[x[1] for x in tissue_sample_pairs]),
        expand("objects/{tissue}_combined.rds", tissue=tissue_set),
        expand("data/pseudobulk/{tissue}_pseudobulk.rds", tissue=tissue_set)
        
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
            -f {input.fragments} \
            -c {input.barcodes} \
            -o {output} \
            -b {input.regions} #\
            #--pic   # there was no --pic option when pseudobulk was first created
        """

rule build_object:
    input:
        done=lambda wc: f"data/{samples[wc.sample]['sample_type']}/.download_complete",
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
        nucleosome_signal_below=lambda wc: samples[wc.sample]["nucleosome_signal"]
    script:
        "code/build_object.R"

rule combine_object:
    input:
        done=lambda wc: [
            f"data/{samples[s]['sample_type']}/.download_complete"
            for s in [k for k, v in samples.items() if v["sample_type"] == wc.tissue]
        ],
        objects=lambda wc: [
            f"objects/{s}_peaks.rds"
            for s in [k for k, v in samples.items() if v["sample_type"] == wc.tissue]
        ]
    output:
        combined="objects/{tissue}_combined.rds",
        integrated="objects/{tissue}_integrated.rds"
    params:
        tissue_name=lambda wc: wc.tissue,
        clustering_res=lambda wc: tissue_clustering_resolution[wc.tissue]
    script:
        "code/combine_objects.R"

# rule annotate_celltypes:
# todo: annotate celltypes for multiome datasets

rule tissue_pseudobulk:
    input:
        tissue_integrated="objects/{tissue}_integrated.rds"
    output:
        pseudobulk="data/pseudobulk/{tissue}_pseudobulk.rds"
    script:
        "code/tissue_pseudobulk.R"

#rule combine_tissue_pseudobulk:
#    input:
#        adrenal="data/pseudobulk/adrenal_pseudobulk.rds",
#        esophagus="data/pseudobulk/esophagus_pseudobulk.rds",
#        heartRV="data/pseudobulk/heartRV_pseudobulk.rds",
#        heart_fetal="data/pseudobulk/heart_fetal_pseudobulk.rds",
#        left_colon="data/pseudobulk/left_colon_pseudobulk.rds",
#        liver="data/pseudobulk/liver_pseudobulk.rds",
#        psoas_muscle="data/pseudobulk/psoas_muscle_pseudobulk.rds"
#    output:
#        tissues_pseudobulk="data/pseudobulk/tissues_pseudobulk.rds"
#    script:
#        "code/combine_pseudobulk.R"

# rule split_matrix


