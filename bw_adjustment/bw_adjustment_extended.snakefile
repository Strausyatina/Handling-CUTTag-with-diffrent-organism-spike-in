# original author of the code: Tarcísio Fontenele de Brito
# changes introduced by Asia Mendelevich
# Extends spikein_norm.chrnames_files.snakefile:
#   - host and spike-in focus chromosomes are both given as chromosome lists
#   - the BAM is reduced to each organism's focus chromosomes, so raw/CPM tracks
#     are computed over that organism's chromosomes ONLY (CPM denominator too)
#   - a blacklist can be given per organism ("none" to skip)
#   - five bigWigs per sample:
#       {sample}.{host}_raw.bw       {sample}.{host}_CPM.bw
#       {sample}.{spikein}_raw.bw    {sample}.{spikein}_CPM.bw
#       {sample}.{host}_scaled.bw    (host, scaled by the spike-in size factor)
#
# Run:
#   bash bw_adjustment_extended.sh ${dir_snakepipes_dnamapping_out} ${path_host_chrs} ${path_spikein_chrs} ${host_name} ${spikein_name} ${blacklist_host} ${blacklist_spikein} ${out_suffix} -j 16

import glob, os

configfile: os.path.join("DNAmapping_organism.yaml") # config (set up locally, relative to defined directory; most probably you have it like that in the snakepipes DNAmapping output directory)

out_suffix = config["out_suffix"]
subdirectory = f"bamCoverage_scaleFactor{out_suffix}"

host = config["host_name"]
spikein = config["spikein_name"]

# "--blackListFileName <bed>", or nothing when the organism has no blacklist.
# An empty value, "none"/"None", or leaving the argument out altogether all mean
# "no blacklist". Note snakemake turns an empty --config value into None.
bl_host = config.get("blacklist_host", "")
bl_spikein = config.get("blacklist_spikein", "")
blacklist_host = "" if not bl_host or str(bl_host).lower() == "none" \
                 else "--blackListFileName " + str(bl_host)
blacklist_spikein = "" if not bl_spikein or str(bl_spikein).lower() == "none" \
                    else "--blackListFileName " + str(bl_spikein)

samplenames, = glob_wildcards(os.path.join("filtered_bam",'{sample}.filtered.bam'))

rule all:
    input:
        os.path.join(subdirectory,"multiBamSummary.spike_in.scaleFactors.tsv"),
        expand(os.path.join(subdirectory, "{sample}." + host + "_raw.bw"), sample=samplenames),
        expand(os.path.join(subdirectory, "{sample}." + host + "_CPM.bw"), sample=samplenames),
        expand(os.path.join(subdirectory, "{sample}." + spikein + "_raw.bw"), sample=samplenames),
        expand(os.path.join(subdirectory, "{sample}." + spikein + "_CPM.bw"), sample=samplenames),
        expand(os.path.join(subdirectory, "{sample}." + host + "_scaled.bw"), sample=samplenames),
        os.path.join(subdirectory, "chromosome_counts.tsv"),
        os.path.join(subdirectory, "multiBamSummary.spike_in.scaleFactors.relative_to_geomean.tsv")

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

# Computing scale factors from spike-in proportions (chrs provided in file):

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

# Size factors for CSAW / DESeq2-style tools:
#
# IMPORTANT -- orientation. multiBamSummary --scalingFactors returns the
# RECIPROCAL of the DESeq2 size factor, because bamCoverage MULTIPLIES coverage
# by it ("The inverse of that is returned, as it's then compatible with
# bamCoverage", deeptools/countReadsPerBin.py::estimateSizeFactors).
# CSAW/DESeq2/edgeR instead DIVIDE counts by the size factor, so the numbers
# have to be inverted before they are used there. Feeding the bamCoverage
# numbers to CSAW unchanged corrects every sample the wrong way round.
#
# Factors are centred on their geometric mean (the edgeR convention), so the
# typical sample sits at 1.0. Centring cancels out of differential testing --
# it only makes the table readable.

rule csaw_size_factors:
    input:
        factortable=rules.scalingFactors.output.factortable,
        counts=os.path.join(subdirectory, "chromosome_counts.tsv"),
        spikein_names=config['path_spikein_chrs']
    output:
        tsv=os.path.join(subdirectory, "multiBamSummary.spike_in.scaleFactors.relative_to_geomean.tsv")
    shell:
        r"""
        # per-sample spike-in read counts, from the idxstats table
        awk -F'\t' 'NR==FNR {{spike[$1]=1; next}} FNR>1 && ($2 in spike) {{reads[$1]+=$4}} \
                    END {{for (s in reads) print s"\t"reads[s]}}' \
            {input.spikein_names} {input.counts} > {output.tsv}.reads.tmp

        # invert the deeptools factors, then centre on the geometric mean
        awk -F'\t' 'NR==FNR {{reads[$1]=$2; next}}
                    FNR==1 {{next}}
                    {{
                      s=$1; sub(/\.filtered\.bam$/, "", s)
                      n++; name[n]=s; dt[n]=$2; size[n]=1/$2; logsum+=log(1/$2)
                    }}
                    END {{
                      g=exp(logsum/n)
                      print "sample\tspikein_reads\tscale_factor_deeptools\tsize_factor\tsize_factor_rel_geomean"
                      for (i=1; i<=n; i++)
                        printf "%s\t%d\t%.6f\t%.6f\t%.6f\n", \
                               name[i], reads[name[i]], dt[i], size[i], size[i]/g
                    }}' \
            {output.tsv}.reads.tmp {input.factortable} > {output.tsv}

        rm -f {output.tsv}.reads.tmp
        """

# Create host-only and spikein-only BAMs (temp: removed once their tracks are built):

rule host_only_bam:
    input:
        bam=os.path.join("filtered_bam", "{sample}.filtered.bam"),
        bai=os.path.join("filtered_bam", "{sample}.filtered.bam.bai"),
        chr_names=config["path_host_chrs"]
    output:
        bam=temp(os.path.join(subdirectory, "host_only_bam", "{sample}.host_only.bam")),
        bai=temp(os.path.join(subdirectory, "host_only_bam", "{sample}.host_only.bam.bai"))
    log:
        os.path.join(subdirectory, "log", "{sample}.host_only_bam.log")
    shell:
        r"""
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        keep=$(samtools view -H {input.bam} | awk -F'\t' '$1=="@SQ"{{for(i=2;i<=NF;i++) if($i~/^SN:/) print substr($i,4)}}' | grep -Fxf {input.chr_names} | tr '\n' ' ')
        samtools view -b {input.bam} $keep > {output.bam} 2> {log}
        samtools index {output.bam} {output.bai} 2>> {log}
        """

rule spikein_only_bam:
    input:
        bam=os.path.join("filtered_bam", "{sample}.filtered.bam"),
        bai=os.path.join("filtered_bam", "{sample}.filtered.bam.bai"),
        chr_names=config["path_spikein_chrs"]
    output:
        bam=temp(os.path.join(subdirectory, "spikein_only_bam", "{sample}.spikein_only.bam")),
        bai=temp(os.path.join(subdirectory, "spikein_only_bam", "{sample}.spikein_only.bam.bai"))
    log:
        os.path.join(subdirectory, "log", "{sample}.spikein_only_bam.log")
    shell:
        r"""
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        keep=$(samtools view -H {input.bam} | awk -F'\t' '$1=="@SQ"{{for(i=2;i<=NF;i++) if($i~/^SN:/) print substr($i,4)}}' | grep -Fxf {input.chr_names} | tr '\n' ' ')
        samtools view -b {input.bam} $keep > {output.bam} 2> {log}
        samtools index {output.bam} {output.bai} 2>> {log}
        """

# Scaling host bw files:

rule bamCoverage_scaleFactor:
    input:
        bam=os.path.join(subdirectory, "host_only_bam", "{sample}.host_only.bam"),
        bai=os.path.join(subdirectory, "host_only_bam", "{sample}.host_only.bam.bai"),
        factortable=rules.scalingFactors.output.factortable
    output:
        bigwig=os.path.join(subdirectory, "{sample}." + host + "_scaled.bw")
    params:
        param_string="--binSize 25 --minMappingQuality 3 --extendReads",
        blacklist=blacklist_host,
        label="{sample}.filtered.bam"   # column 1 of the scale factor table
    threads: 8
    log:
        os.path.join(subdirectory, "log", "bamCoverage.{sample}.scaled.log")
    shell:
        r"""
        bamCoverage -p {threads} {params.param_string} {params.blacklist} \
          -b {input.bam} \
          -o {output.bigwig} \
          --scaleFactor $(awk -F'\t' '$1=="{params.label}" {{print $2}}' {input.factortable}) \
          &> {log}
        """

# Raw and CPM (within host) scaling of bw:

rule bamCoverage_host_raw:
    input:
        bam=os.path.join(subdirectory, "host_only_bam", "{sample}.host_only.bam"),
        bai=os.path.join(subdirectory, "host_only_bam", "{sample}.host_only.bam.bai")
    output:
        bigwig=os.path.join(subdirectory, "{sample}." + host + "_raw.bw"),
        bigwigCPM=os.path.join(subdirectory, "{sample}." + host + "_CPM.bw")
    params:
        param_string="--binSize 25 --minMappingQuality 3 --extendReads",
        blacklist=blacklist_host
    threads: 8
    log:
        os.path.join(subdirectory, "log", "bamCoverage.{sample}.host_raw_cpm.log")
    shell:
        r"""
        bamCoverage -p {threads} {params.param_string} {params.blacklist} \
          -b {input.bam} \
          -o {output.bigwig} \
          &> {log}
        bamCoverage -p {threads} {params.param_string} {params.blacklist} \
          -b {input.bam} \
          -o {output.bigwigCPM} \
          --normalizeUsing CPM \
          &>> {log}
        """

# Raw and CPM (within spike-in) scaling of bw:

rule bamCoverage_spikein_raw:
    input:
        bam=os.path.join(subdirectory, "spikein_only_bam", "{sample}.spikein_only.bam"),
        bai=os.path.join(subdirectory, "spikein_only_bam", "{sample}.spikein_only.bam.bai")
    output:
        bigwig=os.path.join(subdirectory, "{sample}." + spikein + "_raw.bw"),
        bigwigCPM=os.path.join(subdirectory, "{sample}." + spikein + "_CPM.bw")
    params:
        param_string="--binSize 25 --minMappingQuality 3 --extendReads",
        blacklist=blacklist_spikein
    threads: 8
    log:
        os.path.join(subdirectory, "log", "bamCoverage.{sample}.spikein_raw_cpm.log")
    shell:
        r"""
        bamCoverage -p {threads} {params.param_string} {params.blacklist} \
          -b {input.bam} \
          -o {output.bigwig} \
          &> {log}
        bamCoverage -p {threads} {params.param_string} {params.blacklist} \
          -b {input.bam} \
          -o {output.bigwigCPM} \
          --normalizeUsing CPM \
          &>> {log}
        """
