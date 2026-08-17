module load deeptools slurm snakemake bedtools samtools

workdir="$1"
path_host_chrs="$2"
path_spikein_chrs="$3"
out_suffix="$4"
shift 4

snakemake -s spikein_norm.chrnames_files.snakefile -d "$workdir" \
  --latency-wait 120 --keep-going -j 16 \
  "$@" \
  --config \
    path_host_chrs="$path_host_chrs" \
    path_spikein_chrs="$path_spikein_chrs" \
    out_suffix="$out_suffix"

