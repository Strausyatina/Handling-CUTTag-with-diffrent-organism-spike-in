module load deeptools slurm snakemake bedtools samtools

workdir="$1"
path_host_chrs="$2"
path_spikein_chrs="$3"
host_name="$4"
spikein_name="$5"
blacklist_host="$6"
blacklist_spikein="$7"
out_suffix="$8"
shift 8

snakemake -s bw_adjustment_extended.snakefile -d "$workdir" \
  --latency-wait 120 --keep-going -j 16 \
  "$@" \
  --config \
    path_host_chrs="$path_host_chrs" \
    path_spikein_chrs="$path_spikein_chrs" \
    host_name="$host_name" \
    spikein_name="$spikein_name" \
    blacklist_host="$blacklist_host" \
    blacklist_spikein="$blacklist_spikein" \
    out_suffix="$out_suffix"
