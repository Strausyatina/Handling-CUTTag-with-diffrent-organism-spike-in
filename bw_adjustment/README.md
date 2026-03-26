## Rationale

Capturing systemic shifts reflected in spike-in proportions.

## Versioning:

v260325.0 -- almost unchanged, automatically considers `*_spikein` chromosomes as spike-in together, the rest is considered as host. 

`bash spikein_norm.sh ${dir_snakepipes_dnamapping_out} -j 16`

```
example.spikein_norm.sh
> spikein_norm.sh
> spikein_norm.snakefile
```

v260325.1 -- now host and spike-in chromosomes are specified in text files with chromosomes lists (not `_spikein` as before); 
there is also a mandatory suffix for the bamCoverage_scaleFactor output directory now (like `_lambda`)

`bash spikein_norm.chrnames_files.sh ${dir_snakepipes_dnamapping_out} ${path_host_chrs} ${path_spikein_chrs} ${out_suffix} -j 16`

```
example.spikein_norm.chrnames_files.sh
> spikein_norm.chrnames_files.sh
> spikein_norm.chrnames_files.snakefile
```

