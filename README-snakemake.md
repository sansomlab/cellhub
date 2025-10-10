# CellHub Reimplementation

## Reimplementing Annotation in Snakemake

### Prepare environment

```bash
module load snakemake
snakemake --version # 8.4.2
```

### Usage

To generate a `YAML` template:

```bash
snakemake -s ${PATH_TO_SNAKEFILE} --config target=${MODULE_TO_RUN} mode=genconfig --cores 1
```

To make real execution:

```bash
snakemake -s ${PATH_TO_SNAKEFILE} --config target=${MODULE_TO_RUN} --cores ${NCORES} -p --executor ${EXECUTOR} --jobs ${NJOBS}
    # ${EXECUTOR} can be `slurm`
```
