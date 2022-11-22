configfile: "config.yaml"

TREATMENTS = ["addition", "crop", "fert"]
FLUX_PLOT_TYPES = ["", "_cumulative"]

rule targets:
    input:
        expand("figures/{treatment}_plot{plot_type}.png",
            treatment = TREATMENTS, plot_type = FLUX_PLOT_TYPES
        ),
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
        expand("figures/{treatment}_plot{plot_type}.png",
            treatment = TREATMENTS, plot_type = FLUX_PLOT_TYPES
        )
    log:
        log = "logs/generate_flux_plots.txt"
    conda:
        "environment.yml"
    shell:
        """
        {input.r_script} 2> {log}
        """
