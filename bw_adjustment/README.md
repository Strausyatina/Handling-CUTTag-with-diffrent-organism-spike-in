## Rationale

Capturing systemic shifts reflected in spike-in proportions.

## Versioning:


**v260807.0** -- **new `bw_adjustment_extended`**, host, spike-in and spike-in corrected host in one run
(the earlier versions cover the host only). 

`bash bw_adjustment_extended.sh ${dir_snakepipes_dnamapping_out} ${path_host_chrs} ${path_spikein_chrs} ${host_name} ${spikein_name} ${blacklist_host} ${blacklist_spikein} ${out_suffix} -j 16`

```
example.bw_adjustment_extended.sh
> bw_adjustment_extended.sh
> bw_adjustment_extended.snakefile
```

### Details:

Extends `spikein_norm.chrnames_files` to cover host **and** spike-in in one run.
Host and spike-in focus chromosomes are both given as chromosome lists, the BAM
is reduced to each organism's chromosomes, and a blacklist can be given per
organism. A blacklist argument may be an **empty string**, `none`/`None`, or
left out of `--config` altogether -- all mean "no blacklist for this organism",
so the usual case of a host blacklist and no spike-in blacklist is just `""`.
Five bigWigs per sample:

```
{sample}.{host}_raw.bw       {sample}.{host}_CPM.bw
{sample}.{spikein}_raw.bw    {sample}.{spikein}_CPM.bw
{sample}.{host}_scaled.bw    (host, scaled by the spike-in size factor)
```

Track is built from that organism's reduced BAM, so "CPM" counts only
that organism's reads, and the CPM denominator is restricted, not just the output.

#### Size factors for CSAW

`multiBamSummary.spike_in.scaleFactors.relative_to_geomean.tsv`:

| column | meaning |
|---|---|
| `sample` | sample name (no `.filtered.bam` suffix) |
| `spikein_reads` | spike-in reads, summed from `chromosome_counts.tsv` |
| `scale_factor_deeptools` | as produced by `multiBamSummary`; what bamCoverage **multiplies** by |
| `size_factor` | `1 / scale_factor_deeptools`; what CSAW/DESeq2 **divide** by |
| `size_factor_rel_geomean` | `size_factor` centred on its geometric mean (edgeR convention, factors multiply to 1) |


**v260325.2** -- bugfix only, no change to inputs, outputs or file formats.

**v260325.1** -- now host and spike-in chromosomes are specified in text files with chromosomes lists (not `_spikein` as before); 
there is also a mandatory suffix for the bamCoverage_scaleFactor output directory now (like `_lambda`); it will also create a per-chromosome count table for all samples.


`bash spikein_norm.chrnames_files.sh ${dir_snakepipes_dnamapping_out} ${path_host_chrs} ${path_spikein_chrs} ${out_suffix} -j 16`

```
example.spikein_norm.chrnames_files.sh
> spikein_norm.chrnames_files.sh
> spikein_norm.chrnames_files.snakefile
```

**v260325.0** -- almost unchanged, automatically considers `*_spikein` chromosomes as spike-in together, the rest is considered as host. 

`bash spikein_norm.sh ${dir_snakepipes_dnamapping_out} -j 16`

```
example.spikein_norm.sh
> spikein_norm.sh
> spikein_norm.snakefile
```


