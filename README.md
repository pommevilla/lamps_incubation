# LAMPS Incubation

Project repo for analysis done on the crop priming data for the LAMPS incubation project.

## File descriptions

* `results`
    * `stats`
        * `n_anova_results.long.csv` and `n_anova_results.wide.csv`: ANOVA results for NH4, NO3, CO2, and N2O (daily and cumulative), as well as  mineralization and nitrification daily rates (absolute and relative), against Day, Treatment, Addition, and N Rate.
* `figures`
    * `n_figures`
        * [`n_line_plots.png`](figures/n_figures/n_line_plots.png) - daily GHG and N by treatment.
        * [`n_line_plots_cumulative.png`](figures/n_figures/n_line_plots_cumulative.png) - cumulative GHG and N by treatment.
        * [`net_min_nitr_abs_rel.png`](figures/n_figures/net_min_nitr_abs_rel.png) - nitrification and mineralization rates (relative to previous sample and relative to initial sample) by treatment.

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
