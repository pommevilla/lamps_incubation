# LAMPS Incubation

Project repo for analysis done on the crop priming data for the LAMPS incubation project.

## Repeating the analysis

After cloning the repo, restore the conda environment by doing:

```
mamba env create -f environment.yml
```

Then, activate the environment:

```
conda activate lamps_incubation
```

You can then run the pipeline with:

```
snakemake -c 1
```

_(You can also do `snakemake -c 2` or however many cores you'd like to use.)_

This will restore the R environment, regenerate all the figures, and render the website.
