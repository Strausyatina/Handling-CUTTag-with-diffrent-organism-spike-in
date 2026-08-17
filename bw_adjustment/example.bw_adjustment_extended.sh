# Arguments:
bash bw_adjustment_extended.sh ${dir_snakepipes_dnamapping_out} ${path_host_chrs} ${path_spikein_chrs} ${host_name} ${spikein_name} ${blacklist_host} ${blacklist_spikein} ${out_suffix} -j 16
#
#   1. dir_snakepipes_dnamapping_out  snakePipes DNAmapping output directory; becomes
#                                     snakemake -d, so all output is written inside it.
#                                     Must contain filtered_bam/<sample>.filtered.bam
#                                     (+ .bam.bai) -- sample names are taken from those
#                                     filenames -- and DNAmapping_organism.yaml, which is
#                                     read for genome_index (the .fai).
#   2. path_host_chrs                 text file, one host chromosome name per line
#   3. path_spikein_chrs              text file, one spike-in chromosome name per line
#                                     (names as they appear in the hybrid BAM, e.g.
#                                     2L_spikein). Also defines the regions the spike-in
#                                     size factors are computed over.
#                                     For 2. and 3.: names must match the BAM header
#                                     exactly; names absent from the BAM are ignored, so
#                                     one list can serve several references.
#                                     These are FOCUS chromosomes -- the ones to analyse,
#                                     not every chromosome of that organism. Typically the
#                                     autosomes (chr1..chr22), or autosomes + chrX. What is
#                                     left out is left out on purpose: sex chromosomes when
#                                     the cohort mixes sexes or copy number differs, chrM
#                                     (tiny, extremely high coverage, distorts CPM), and
#                                     unplaced/random/alt contigs (unreliable mapping).
#                                     The choice matters twice over: it decides which
#                                     regions end up in the track AND which reads count
#                                     towards the CPM denominator, since the BAM is reduced
#                                     to these chromosomes before bamCoverage runs.
#                                     Supplied lists live in support_files/.
#   4. host_name                      label used in the output filenames, e.g. human.
#   5. spikein_name                   label used in the output filenames, e.g. S2.
#                                     4. and 5. are labels only -- nothing checks that
#                                     they agree with the chromosome lists.
#   6. blacklist_host                 BED of regions to exclude from the host tracks;
#                                     "" or none for no blacklist.
#   7. blacklist_spikein              same for the spike-in tracks; usually "" or none,
#                                     as spike-in genomes rarely have a blacklist.
#   8. out_suffix                     suffix for the output directory
#                                     bamCoverage_scaleFactor${out_suffix}. Give each run
#                                     its own, so two runs do not overwrite one another.
#                                     A leading underscore is conventional, not automatic.
#   [snakemake args...]               passed straight to snakemake, and placed after this
#                                     script's defaults so they win: -j N (cores, default
#                                     16), -n (dry run), -p (print commands),
#                                     --rerun-incomplete, --unlock.
#
# One chromosome name per line -- this is required, not a convention.
#
# Chromosome-name extensions must be consistent everywhere:
#   In a hybrid genome the spike-in chromosomes usually carry an extension, e.g.
#   2L -> 2L_spikein (snakePipes writes this with --spikeinExt, default _spikein).
#   Whatever extension is used, the SAME spelling has to appear in:
#     - the BAM header            (@SQ SN:2L_spikein)
#     - the chromosome lists      (2. path_host_chrs, 3. path_spikein_chrs)
#     - the blacklist BED files   (6. blacklist_host, 7. blacklist_spikein, column 1)
#     - the genome_index (.fai) named by DNAmapping_organism.yaml, which is what
#       the spike-in regions for the size factors are built from.
#
#   If a blacklist seems to do nothing, check its naming first:
#         cut -f1 <blacklist.bed> | sort -u
#         samtools view -H <bam> | awk -F'\t' '$1=="@SQ"{print $2}' | sort -u
#
# Output, into ${dir_snakepipes_dnamapping_out}/bamCoverage_scaleFactor${out_suffix}/ :
#   <sample>.${host_name}_raw.bw       raw,  host chrs only
#   <sample>.${host_name}_CPM.bw       CPM,  host chrs only (host reads only in the denominator)
#   <sample>.${spikein_name}_raw.bw    raw,  spike-in chrs only
#   <sample>.${spikein_name}_CPM.bw    CPM,  spike-in chrs only
#   <sample>.${host_name}_scaled.bw    host, scaled by the spike-in size factor
#   multiBamSummary.spike_in.scaleFactors.tsv                    factors for bamCoverage
#   multiBamSummary.spike_in.scaleFactors.relative_to_geomean.tsv   factors for CSAW
#   chromosome_counts.tsv                                        per-chromosome read counts

# Example: human host + S2 (dm6) spike-in, blacklist for the host only.
# A blacklist argument may be "" or "none" - both mean no blacklist:
# bash bw_adjustment_extended.sh ${dir_snakepipes_dnamapping_out} \
#   /support_files/GRCh38_chromosomes.focus.txt \
#   /support_files/dm6_chromosomes.focus_spikein.txt \
#   human S2 \
#   /path/to/human_blacklist.bed none \
#   _human_s2 -j 16
