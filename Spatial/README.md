# SpaISO-Seq pipeline

## Prepare

- [sicelore](https://github.com/ucagenomix/sicelore)
- [minimap2](https://github.com/lh3/minimap2)
- [poa](https://github.com/tanghaibao/bio-pipeline/tree/master/poaV2)
- samtools

<font color=red>Attention</font>

Rename poa to spoa and make sure `spoa` in your environment variable.

## Run

```shell
# Dry run
snakemake -s sicelorel_pipe.smk --np
```

