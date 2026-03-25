# original author of the code: Yinxiu Zhan  
# (minor) changes introduced by Asia Mendelevich 
# Changes: 
#   - multiBamSummary bins --> multiBamSummary BED-file (`bins --region` accepted only 1 region, we may have several)
#
# Run:
#   spikein_norm.sh:
#     module load deeptools slurm snakemake bedtools
#     read workdir snakemake_options <<< "$@"
#     snakemake -s spikein_norm.snakefile -d $workdir $snakemake_options --latency-wait 120 --keep-going --cores 16 
#   
#   bash spikein_norm.sh ${dir_snakepipes_dnamapping_out} -j 16

import glob, os

configfile: os.path.join("DNAmapping_organism.yaml") # config (set up locally, relative to defined directory; most probably you have it like that in the snakepipes DNAmapping output directory)

samplenames, = glob_wildcards(os.path.join("filtered_bam",'{sample}.filtered.bam'))

rule all:
    input:
        os.path.join("bamCoverage_scaleFactor","multiBamSummary.spike_in.scaleFactors.tsv"),
        expand(os.path.join("bamCoverage_scaleFactor","{sample}.scaled.bw"), sample = samplenames)

# Computing scale factors from spike-in proportions (defined: *_spikein chrs):

rule spike_in_region:
    input:
        fai=config['genome_index']
    output:
        bed=os.path.join("bamCoverage_scaleFactor","spike_in.bed"),
        bins=os.path.join("bamCoverage_scaleFactor","spike_in.bins.bed")
    params:
        binsize=1000
    shell:
        r"""
        grep '_spikein' {input.fai} | awk -v OFS='\t' '{{print($1,0,$2)}}' > {output.bed}; 
        bedtools makewindows -b {output.bed} -w {params.binsize} > {output.bins}
        """

rule scalingFactors:
    input:
        bam_files=expand(os.path.join("filtered_bam","{sample}.filtered.bam"), sample = samplenames),
        bai_files=expand(os.path.join("filtered_bam","{sample}.filtered.bam.bai"), sample = samplenames),
        spikein_bed=rules.spike_in_region.output.bins
    output:
        outtable=os.path.join("bamCoverage_scaleFactor","multiBamSummary.spike_in.npz"),
        rawtable=os.path.join("bamCoverage_scaleFactor","multiBamSummary.spike_in.rawcounts.tsv"),
        factortable=os.path.join("bamCoverage_scaleFactor","multiBamSummary.spike_in.scaleFactors.tsv")
    threads: 32
    log: os.path.join("bamCoverage_scaleFactor", "log", "multiBamSummary.bedfile.log")
    shell:
        "multiBamSummary BED-file --BED {input.spikein_bed} -p {threads} --bamfiles {input.bam_files} --outFileName {output.outtable} --scalingFactors {output.factortable} --outRawCounts {output.rawtable} &> {log}"

# Scaling host bw files:

rule bamCoverage_scaleFactor:
    input:
        bam=os.path.join("filtered_bam","{sample}.filtered.bam"),
        factortable=rules.scalingFactors.output.factortable
    output:
        bigwig=os.path.join("bamCoverage_scaleFactor", "{sample}.scaled.bw")
    params:
        param_string="--binSize 25 --minMappingQuality 3 --extendReads",
        pattern="{sample}.filtered.bam"
    threads: 8
    log: os.path.join("bamCoverage_scaleFactor", "log", "bamCoverage.{sample}.scaled.log")
    shell:
        "bamCoverage -p {threads} {params.param_string} -b {input.bam} -o {output.bigwig} --scaleFactor $(grep '{params.pattern}' {input.factortable} | cut -f 2) &> {log}"
