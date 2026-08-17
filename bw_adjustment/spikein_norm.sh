module load deeptools slurm snakemake bedtools

read workdir snakemake_options <<< "$@"

snakemake -s spikein_norm.snakefile -d $workdir --latency-wait 120 --keep-going -j 16 $snakemake_options
