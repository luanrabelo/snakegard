

def get_multiqc_files_by_species(wildcards):
    """Gathers all QC files for a specific species for its MultiQC report."""
    qc_files = []
    species_wc = wildcards.species
    for entry in config["samples"]:
        species = entry["species"].replace(" ", "_")
        if species == species_wc:
            for sample in entry["sra"]:
                sra_id = sample["id"]
                qc_files.extend([
                    f"results/02-fastqc_raw/{species}/{sra_id}_R1_fastqc.zip",
                    f"results/02-fastqc_raw/{species}/{sra_id}_R2_fastqc.zip",
                    f"results/03-trim_data/{species}/fastp/{sra_id}.fastp.html"
                    f"results/03-trim_data/{species}/fastp/{sra_id}.fastp.json",
                    f"results/04-fastqc_trimmed/{species}/{sra_id}_R1.trimmed_fastqc.zip",
                    f"results/04-fastqc_trimmed/{species}/{sra_id}_R2.trimmed_fastqc.zip"
                ])
    
    return qc_files

rule multiqc_by_species:
    """Aggregates all QC results for a given species into a single report."""
    input:
        get_multiqc_files_by_species(wildcards)
    output:
        report="results/05-multiqc/{species}/multiqc_report.html",
        data_dir="results/05-multiqc/{species}/multiqc_data"
    params:
        analysis_dirs=lambda wildcards: [
            f"results/02-fastqc_raw/{wildcards.species}",
            f"results/03-trim_data/logs/{wildcards.species}",
            f"results/04-fastqc_trimmed/{wildcards.species}"
        ],
        outdir="results/05-multiqc/{species}"
    threads: 1
    run:
        logger = setup_logger(f"multiqc_{wildcards.species}", log_file)
        try:
            logger.info(f"Aggregating QC reports for species {wildcards.species}.")
            cmd = ["multiqc", "--force", "--outdir", params.outdir] + params.analysis_dirs
            subprocess.run(cmd, check=True)
            logger.info(f"MultiQC report for {wildcards.species} generated successfully in {params.outdir}")
        except subprocess.CalledProcessError as e:
            logger.error(f"MultiQC failed for {wildcards.species}: {e}")
            raise e