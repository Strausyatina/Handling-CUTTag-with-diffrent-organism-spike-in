module load deeptools slurm snakemake bedtools

read workdir snakemake_options <<< "$@"

snakemake -s spikein_norm.snakefile -d $workdir $snakemake_options --latency-wait 120 --keep-going --cores 16
