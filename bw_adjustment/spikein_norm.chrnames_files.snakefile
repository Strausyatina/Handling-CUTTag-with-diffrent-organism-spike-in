# original author of the code: Yinxiu Zhan  
# changes introduced by Asia Mendelevich 
# Changes: 
#   - multiBamSummary bins --> multiBamSummary BED-file (`bins --region` accepted only 1 region, we may have several)
#   - now host and spike-in chromosomes are specified in text files with chromosomes lists (not `_spikein` as before) 
#   - there is also a mandatory suffix for the bamCoverage_scaleFactor output directory now (like `_lambda`)
#   - it will also create a per-chromosome raw count table
#
# Run:
#   bash spikein_norm.chrnames_files.sh ${dir_snakepipes_dnamapping_out} ${path_host_chrs} ${path_spikein_chrs} ${out_suffix} -j 16

import glob, os

configfile: os.path.join("DNAmapping_organism.yaml") # config (set up locally, relative to defined directory; most probably you have it like that in the snakepipes DNAmapping output directory)

out_suffix = config["out_suffix"]
subdirectory = f"bamCoverage_scaleFactor{out_suffix}"

samplenames, = glob_wildcards(os.path.join("filtered_bam",'{sample}.filtered.bam'))

rule all:
    input:
        os.path.join(subdirectory,"multiBamSummary.spike_in.scaleFactors.tsv"),
        expand(os.path.join(subdirectory,"{sample}.scaled.bw"), sample = samplenames),
        os.path.join(subdirectory, "chromosome_counts.tsv")

# Compute per-chromosome counts table:

rule chromosome_counts_per_sample:
    input:
        bam=os.path.join("filtered_bam", "{sample}.filtered.bam"),
        bai=os.path.join("filtered_bam", "{sample}.filtered.bam.bai")
    output:
        tsv=os.path.join(subdirectory, "chrom_counts", "{sample}.idxstats.tsv")
    log:
        os.path.join(subdirectory, "log", "{sample}.idxstats.log")
    shell:
        r"""
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        samtools idxstats {input.bam} \
          | awk 'BEGIN{{OFS="\t"; print "sample","chromosome","length","mapped","unmapped"}} {{print "{wildcards.sample}", $1, $2, $3, $4}}' \
          > {output.tsv} 2> {log}
        """

rule chromosome_counts_table:
    input:
        expand(os.path.join(subdirectory, "chrom_counts", "{sample}.idxstats.tsv"), sample=samplenames)
    output:
        os.path.join(subdirectory, "chromosome_counts.tsv")
    shell:
        r"""
        awk 'FNR==1 && NR!=1 {{next}} {{print}}' {input} > {output}
        """

# Computing scale factors from spike-in proportions (defined: *_spikein chrs):

rule spike_in_region:
    input:
        fai=config['genome_index'],
        spikein_names=config['path_spikein_chrs']
    output:
        bed=os.path.join(subdirectory,"spike_in.bed"),
        bins=os.path.join(subdirectory,"spike_in.bins.bed")
    params:
        binsize=1000
    shell:
        r"""
        awk 'BEGIN{{OFS="\t"}} NR==FNR {{keep[$1]=1; next}} ($1 in keep) {{print $1, 0, $2}}' {input.spikein_names} {input.fai} > {output.bed};
        bedtools makewindows -b {output.bed} -w {params.binsize} > {output.bins}
        """

rule scalingFactors:
    input:
        bam_files=expand(os.path.join("filtered_bam","{sample}.filtered.bam"), sample = samplenames),
        bai_files=expand(os.path.join("filtered_bam","{sample}.filtered.bam.bai"), sample = samplenames),
        spikein_bed=rules.spike_in_region.output.bins
    output:
        outtable=os.path.join(subdirectory,"multiBamSummary.spike_in.npz"),
        rawtable=os.path.join(subdirectory,"multiBamSummary.spike_in.rawcounts.tsv"),
        factortable=os.path.join(subdirectory,"multiBamSummary.spike_in.scaleFactors.tsv")
    threads: 32
    log: os.path.join(subdirectory, "log", "multiBamSummary.bedfile.log")
    shell:
        "multiBamSummary BED-file --BED {input.spikein_bed} -p {threads} --bamfiles {input.bam_files} --outFileName {output.outtable} --scalingFactors {output.factortable} --outRawCounts {output.rawtable} &> {log}"

# Create host-only BAMs:

rule host_only_bam:
    input:
        bam=os.path.join("filtered_bam", "{sample}.filtered.bam"),
        bai=os.path.join("filtered_bam", "{sample}.filtered.bam.bai"),
        host_names=config["path_host_chrs"]
    output:
        bam=os.path.join(subdirectory, "host_only_bam", "{sample}.host_only.bam"),
        bai=os.path.join(subdirectory, "host_only_bam", "{sample}.host_only.bam.bai")
    log:
        os.path.join(subdirectory, "log", "{sample}.host_only_bam.log")
    shell:
        r"""
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        samtools view -b {input.bam} $(cat {input.host_names}) > {output.bam} 2> {log}
        samtools index {output.bam} {output.bai} 2>> {log}
        """

# Scaling host bw files:

rule bamCoverage_scaleFactor:
    input:
        bam=os.path.join(subdirectory, "host_only_bam", "{sample}.host_only.bam"),
        bai=os.path.join(subdirectory, "host_only_bam", "{sample}.host_only.bam.bai"),
        factortable=rules.scalingFactors.output.factortable
    output:
        bigwig=os.path.join(subdirectory, "{sample}.scaled.bw")
    params:
        param_string="--binSize 25 --minMappingQuality 3 --extendReads",
        pattern="{sample}.filtered.bam"
    threads: 8
    log:
        os.path.join(subdirectory, "log", "bamCoverage.{sample}.scaled.log")
    shell:
        r"""
        bamCoverage -p {threads} {params.param_string} \
          -b {input.bam} \
          -o {output.bigwig} \
          --scaleFactor $(grep '{params.pattern}' {input.factortable} | cut -f 2) \
          &> {log}
        """

