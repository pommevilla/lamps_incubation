rule targets:
    input:
        "figures/addition_plot_cumulative.png",
        "figures/addition_plot.png",
        "figures/crop_plot_cumulative.png",
        "figures/crop_plot.png",
        "figures/fert_plot_cumulative.png",
        "figures/fert_plot.png",
        "figures/snakemake_dag.png"

rule generate_snakemake_dag:
    input:
        script = "code/make_snakemake_dag.sh"
    output:
        "figures/snakemake_dag.png"
    shell:
        """
        {input.script}
        """

rule generate_flux_and_amoa_abundance_plots:
    input:
        r_script = "code/create_line_charts.R"
    output:
        "figures/addition_plot_cumulative.png",
        "figures/addition_plot.png",
        "figures/crop_plot_cumulative.png",
        "figures/crop_plot.png",
        "figures/fert_plot_cumulative.png",
        "figures/fert_plot.png"
    conda:
        "environment.yml"
    shell:
        """
        {input.r_script}
        """
