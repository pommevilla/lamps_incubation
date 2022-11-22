rule targets:
    input:
        "figures/addition_plot_cumulative.png",
        "figures/addition_plot.png",
        "figures/crop_plot_cumulative.png",
        "figures/crop_plot.png",
        "figures/fert_plot_cumulative.png",
        "figures/fert_plot.png",
        "figures/snakemake_dag.png",
        ".lamps_incubation_renv_restored"

rule restore_renv:
    input:
        r_script = "code/restore_renv.R"
    output:
        touch(".lamps_incubation_renv_restored")
    log:
        log = "logs/restore_r.txt"
    conda:
        "environment.yml"
    shell:
        """
        {input.r_script} 2> {log}
        """

rule generate_snakemake_dag:
    input:
        script = "code/make_snakemake_dag.sh"
    output:
        "figures/snakemake_dag.png"
    log:
        log = "logs/generate_snakemake_dag.txt"
    conda:
        "environment.yml"
    shell:
        """
        {input.script} 2> {log}
        """

rule generate_flux_and_amoa_abundance_plots:
    input:
        restored_env = ".lamps_incubation_renv_restored",
        r_script = "code/create_line_charts.R"
    output:
        "figures/addition_plot_cumulative.png",
        "figures/addition_plot.png",
        "figures/crop_plot_cumulative.png",
        "figures/crop_plot.png",
        "figures/fert_plot_cumulative.png",
        "figures/fert_plot.png"
    log:
        log = "logs/generate_flux_plots.txt"
    conda:
        "environment.yml"
    shell:
        """
        {input.r_script} 2> {log}
        """
