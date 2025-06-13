# snakegard/Snakefile

import logging
import os
import subprocess
import tempfile
import shutil
import time
from snakemake.utils import min_version

# ===================================================================
# --- 1. Global Configuration & Helper Functions
# ===================================================================

configfile: "config.yaml"
min_version("7.25.0")

def get_required_containers():
    """Helper function to find which containers are needed based on the config."""
    needed_containers = set()
    for entry in config["samples"]:
        for sample in entry["sra"]:
            assembler = sample.get("assembler")
            if assembler in config["containers"]:
                needed_containers.add(assembler)
    return [os.path.join(config["container_dir"], f"{name}.sif") for name in needed_containers]

def get_all_target_files():
    """
    Collect all final files this pipeline should produce.
    This now includes both the final assembly and the aggregate QC report for each species.
    """
    targets = get_required_containers()
    
    for entry in config["samples"]:
        species = entry["species"].replace(" ", "_")
        
        # Add the per-species QC report as a desired output
        targets.append(f"results/05-multiqc/{species}/multiqc_report.html")
        
        for sample in entry["sra"]:
            sra_id = sample["id"]
            reference = sample.get("reference")
            mapper = sample.get("mapper", "bowtie2")
            assembler = sample.get("assembler", "spades")
            # Corrected, unambiguous output path for final assembly
            targets.append(f"results/07-assembly/{species}/{sra_id}/{reference}/{mapper}/{assembler}/contigs.fasta")
            
    return targets

def get_qc_files_by_species(wildcards):
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
                    f"results/03-trim_data/logs/{species}/{sra_id}.fastp.json",
                    f"results/04-fastqc_trimmed/{species}/{sra_id}_R1.trimmed_fastqc.zip",
                    f"results/04-fastqc_trimmed/{species}/{sra_id}_R2.trimmed_fastqc.zip"
                ])
    return qc_files

def get_mapper_and_reference(wildcards):
    """Helper function to find the correct mapper, reference, and assembler for a given sample."""
    for entry in config["samples"]:
        if entry["species"].replace(" ", "_") == wildcards.species:
            for sample in entry["sra"]:
                if sample["id"] == wildcards.sra:
                    return {
                        "mapper": sample.get("mapper", "bowtie2"),
                        "reference": sample.get("reference"),
                        "assembler": sample.get("assembler", "spades")
                    }
    return {}

def setup_logger(name, log_file):
    """Function to configure a logger to write to a file and the console."""
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    if logger.hasHandlers():
        logger.handlers.clear()
    
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    return logger

# ===================================================================
# --- 2. Main Rule & Pipeline Definition
# ===================================================================

rule all:
    input:
        get_all_target_files()

# --- Rules for Downloading ---
rule download_container:
    """
    Downloads and prepares a Singularity container (.sif file).
    This version is extra safe: it downloads to a temporary path and only
    renames upon successful completion to avoid corrupted files.
    """
    output:
        sif=os.path.join(config["container_dir"], "{container_name}.sif")
    params:
        uri=lambda wildcards: config["containers"][wildcards.container_name]
    threads: 1
    run:
        log_file = str(output.sif) + ".log"
        logger = setup_logger(f"download_container_{wildcards.container_name}", log_file)
        
        # Define final and temporary paths
        final_path = str(output.sif)
        tmp_path = final_path + ".tmp"

        try:
            logger.info(f"Starting job for container '{wildcards.container_name}'.")
            if os.path.exists(final_path):
                logger.info(f"Container already exists at {final_path}. Skipping.")
                return

            os.makedirs(os.path.dirname(final_path), exist_ok=True)
            uri = params.uri
            user_agent = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.36"

            # Download logic now uses tmp_path
            if "dl=1" in uri:
                # ... (mesma lógica do ZIP, mas baixando para tmp_path)
                pass # A lógica completa já está no seu arquivo
            elif uri.startswith("http"):
                logger.info(f"Starting direct download to temporary file: {tmp_path}")
                subprocess.run(["wget", "--quiet", "--user-agent", user_agent, "-O", tmp_path, uri], check=True)
            else:
                logger.info(f"Pulling container to temporary file: {tmp_path}")
                cmd = ["singularity", "pull", "--name", tmp_path, uri]
                # ... (resto da lógica do singularity pull)
            
            # Atomic move from temporary to final path
            logger.info(f"Download complete. Renaming {tmp_path} to {final_path}")
            shutil.move(tmp_path, final_path)

            logger.info(f"Job for container '{wildcards.container_name}' finished successfully.")
        except Exception as e:
            logger.error(f"Job for container '{wildcards.container_name}' failed: {e}")
            # Cleanup temporary file on error
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
            raise e

rule download_reference:
    """Downloads a nucleotide sequence from NCBI using its accession number."""
    output:
        fasta="results/00-references/{reference}.fasta"
    params:
        db="nucleotide"
    threads: 1
    run:
        log_file = str(output.fasta) + ".log"
        logger = setup_logger(f"download_reference_{wildcards.reference}", log_file)
        try:
            logger.info(f"Starting job for reference '{wildcards.reference}'.")
            if os.path.exists(str(output.fasta)):
                logger.info("Reference file already exists. Skipping.")
            else:
                max_retries = 3
                retry_wait_seconds = 20
                downloaded_content = ""
                for attempt in range(max_retries):
                    logger.info(f"Downloading sequence (Attempt {attempt + 1}/{max_retries})...")
                    try:
                        cmd = ["efetch", "-db", params.db, "-id", wildcards.reference, "-format", "fasta"]
                        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=300)
                        if result.stdout and not result.stdout.isspace():
                            downloaded_content = result.stdout
                            logger.info("Successfully downloaded content.")
                            break
                        else:
                            logger.warning(f"Attempt {attempt + 1} resulted in empty content from NCBI.")
                    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
                        logger.warning(f"Attempt {attempt + 1} failed with an error: {e}")
                    if attempt + 1 < max_retries:
                        logger.info(f"Waiting {retry_wait_seconds} seconds before retrying.")
                        time.sleep(retry_wait_seconds)
                if not downloaded_content:
                    raise RuntimeError(f"Downloaded content for {wildcards.reference} is empty after {max_retries} attempts.")
                with open(str(output.fasta), "w") as f:
                    f.write(downloaded_content)
            logger.info(f"Job for reference '{wildcards.reference}' finished successfully.")
        except Exception as e:
            logger.error(f"Job for reference '{wildcards.reference}' failed: {e}")
            if os.path.exists(str(output.fasta)):
                os.remove(str(output.fasta))
            raise e

rule download_sra:
    """Downloads paired-end FASTQ data from SRA for a given accession number."""
    output:
        r1=temp("results/01-raw_data/{species}/{sra}_R1.fastq"),
        r2=temp("results/01-raw_data/{species}/{sra}_R2.fastq")
    threads: 8
    run:
        log_file = f"results/01-raw_data/logs/{wildcards.species}_{wildcards.sra}_download.log"
        logger = setup_logger(f"download_sra_{wildcards.sra}", log_file)
        try:
            logger.info(f"Starting job for SRA accession '{wildcards.sra}'.")
            if os.path.exists(str(output.r1)) and os.path.exists(str(output.r2)):
                 logger.info("Temporary FASTQ files already exist. Skipping download.")
            else:
                with tempfile.TemporaryDirectory() as tmp_dir:
                    logger.info(f"Downloading {wildcards.sra} to temporary directory {tmp_dir} using {threads} threads...")
                    cmd = ["fasterq-dump", wildcards.sra, "--split-files", "--threads", str(threads), "--outdir", tmp_dir]
                    subprocess.run(cmd, check=True)
                    src_r1 = os.path.join(tmp_dir, f"{wildcards.sra}_1.fastq")
                    src_r2 = os.path.join(tmp_dir, f"{wildcards.sra}_2.fastq")
                    if not os.path.exists(src_r1) or not os.path.exists(src_r2):
                        raise FileNotFoundError(f"fasterq-dump did not produce the expected FASTQ files in {tmp_dir}")
                    logger.info("Download complete. Moving files to final destination.")
                    shutil.move(src_r1, str(output.r1))
                    shutil.move(src_r2, str(output.r2))
            logger.info(f"Job for SRA accession '{wildcards.sra}' finished successfully.")
        except Exception as e:
            logger.error(f"Job for SRA accession '{wildcards.sra}' failed: {e}")
            if os.path.exists(str(output.r1)): os.remove(str(output.r1))
            if os.path.exists(str(output.r2)): os.remove(str(output.r2))
            raise e

rule compress_fastq:
    """Compresses the downloaded FASTQ files using gzip or pigz."""
    input:
        r1="results/01-raw_data/{species}/{sra}_R1.fastq",
        r2="results/01-raw_data/{species}/{sra}_R2.fastq"
    output:
        r1_gz="results/01-raw_data/{species}/{sra}_R1.fastq.gz",
        r2_gz="results/01-raw_data/{species}/{sra}_R2.fastq.gz"
    threads: 2
    shell:
        """
        if command -v pigz &> /dev/null; then
            pigz -p {threads} -c {input.r1} > {output.r1_gz}
            pigz -p {threads} -c {input.r2} > {output.r2_gz}
        else
            gzip -c {input.r1} > {output.r1_gz}
            gzip -c {input.r2} > {output.r2_gz}
        fi
        """

# --- Rules for QC and Trimming ---
rule fastqc_raw:
    """Runs FastQC on the raw, compressed FASTQ files."""
    input:
        r1="results/01-raw_data/{species}/{sra}_R1.fastq.gz",
        r2="results/01-raw_data/{species}/{sra}_R2.fastq.gz"
    output:
        r1_zip="results/02-fastqc_raw/{species}/{sra}_R1_fastqc.zip",
        r2_zip="results/02-fastqc_raw/{species}/{sra}_R2_fastqc.zip"
    params:
        outdir="results/02-fastqc_raw/{species}"
    threads: 2
    run:
        log_file = f"{params.outdir}/logs/{wildcards.sra}_fastqc_raw.log"
        logger = setup_logger(f"fastqc_raw_{wildcards.sra}", log_file)
        try:
            logger.info(f"Running FastQC on raw reads for {wildcards.sra}.")
            cmd = ["fastqc", "--threads", str(threads), "--outdir", params.outdir, input.r1, input.r2]
            subprocess.run(cmd, check=True, capture_output=True)
            logger.info(f"FastQC on raw reads for {wildcards.sra} completed successfully.")
        except subprocess.CalledProcessError as e:
            logger.error(f"FastQC on raw reads for {wildcards.sra} failed: {e.stderr.decode()}")
            raise e

rule trim_reads:
    """Uses fastp to perform adapter trimming and quality filtering."""
    input:
        r1="results/01-raw_data/{species}/{sra}_R1.fastq.gz",
        r2="results/01-raw_data/{species}/{sra}_R2.fastq.gz"
    output:
        r1_trim="results/03-trim_data/{species}/{sra}_R1.trimmed.fastq.gz",
        r2_trim="results/03-trim_data/{species}/{sra}_R2.trimmed.fastq.gz",
        json="results/03-trim_data/logs/{species}/{sra}.fastp.json",
        html="results/03-trim_data/logs/{species}/{sra}.fastp.html"
    params:
        extra=config.get("params", {}).get("fastp", "")
    threads: 8
    run:
        log_file = f"results/03-trim_data/logs/{wildcards.species}/{wildcards.sra}_fastp.log"
        logger = setup_logger(f"trim_reads_{wildcards.sra}", log_file)
        try:
            logger.info(f"Running fastp for {wildcards.sra} with {threads} threads.")
            os.makedirs(os.path.dirname(str(output.json)), exist_ok=True)
            cmd = ["fastp", "--in1", input.r1, "--in2", input.r2, "--out1", output.r1_trim, "--out2", output.r2_trim, "--json", output.json, "--html", output.html, "--thread", str(threads)]
            if isinstance(params.extra, dict):
                for k, v in params.extra.items():
                    cmd.append(str(k))
                    if v is not True:
                        cmd.append(str(v))
            elif isinstance(params.extra, str) and params.extra:
                cmd.extend(params.extra.split())
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"fastp for {wildcards.sra} completed successfully.")
        except subprocess.CalledProcessError as e:
            logger.error(f"fastp for {wildcards.sra} failed:\n{e.stderr}")
            raise e

rule fastqc_trimmed:
    """Runs FastQC on the trimmed FASTQ files."""
    input:
        r1="results/03-trim_data/{species}/{sra}_R1.trimmed.fastq.gz",
        r2="results/03-trim_data/{species}/{sra}_R2.trimmed.fastq.gz"
    output:
        r1_zip="results/04-fastqc_trimmed/{species}/{sra}_R1.trimmed_fastqc.zip",
        r2_zip="results/04-fastqc_trimmed/{species}/{sra}_R2.trimmed_fastqc.zip"
    params:
        outdir="results/04-fastqc_trimmed/{species}"
    threads: 2
    run:
        log_file = f"{params.outdir}/logs/{wildcards.sra}_fastqc_trimmed.log"
        logger = setup_logger(f"fastqc_trimmed_{wildcards.sra}", log_file)
        try:
            logger.info(f"Running FastQC on trimmed reads for {wildcards.sra}.")
            cmd = ["fastqc", "--threads", str(threads), "--outdir", params.outdir, input.r1, input.r2]
            subprocess.run(cmd, check=True, capture_output=True)
            logger.info(f"FastQC on trimmed reads for {wildcards.sra} completed successfully.")
        except subprocess.CalledProcessError as e:
            logger.error(f"FastQC on trimmed reads for {wildcards.sra} failed: {e.stderr.decode()}")
            raise e

rule multiqc_by_species:
    """Aggregates all QC results for a given species into a single report."""
    input:
        get_qc_files_by_species
    output:
        report=report("results/05-multiqc/{species}/multiqc_report.html", caption="../report/multiqc.rst", category="Aggregate QC"),
        data_dir=directory("results/05-multiqc/{species}/multiqc_data")
    params:
        analysis_dirs=lambda wildcards: [
            f"results/02-fastqc_raw/{wildcards.species}",
            f"results/03-trim_data/logs/{wildcards.species}",
            f"results/04-fastqc_trimmed/{wildcards.species}"
        ],
        outdir="results/05-multiqc/{species}"
    threads: 1
    run:
        log_file = os.path.join(params.outdir, "multiqc.log")
        logger = setup_logger(f"multiqc_{wildcards.species}", log_file)
        try:
            logger.info(f"Aggregating QC reports for species {wildcards.species}.")
            cmd = ["multiqc", "--force", "--outdir", params.outdir] + params.analysis_dirs
            subprocess.run(cmd, check=True)
            logger.info(f"MultiQC report for {wildcards.species} generated successfully in {params.outdir}")
        except subprocess.CalledProcessError as e:
            logger.error(f"MultiQC failed for {wildcards.species}: {e}")
            raise e

# --- Rules for Mapping and Assembly ---
rule index_reference:
    """Creates an index of the reference sequence for a specific mapper."""
    input:
        ref="results/00-references/{reference}.fasta"
    output:
        done=touch("results/06-mapping/index/{reference}_{mapper}/index.done")
    params:
        ref_prefix="results/06-mapping/index/{reference}_{mapper}/{reference}",
        ref_dir="results/06-mapping/index/{reference}_{mapper}/"
    threads: 8
    run:
        logger = setup_logger(f"index_{wildcards.mapper}_{wildcards.reference}", os.path.join(params.ref_dir, "index.log"))
        try:
            mapper = wildcards.mapper
            logger.info(f"Creating index for {input.ref} using {mapper}.")
            os.makedirs(params.ref_dir, exist_ok=True)
            cmd = []
            if mapper == "bowtie2":
                cmd = ["bowtie2-build", "--threads", str(threads), input.ref, params.ref_prefix]
            elif mapper == "bwa":
                cmd = ["bwa", "index", "-p", params.ref_prefix, input.ref]
            elif mapper == "star":
                cmd = ["STAR", "--runThreadN", str(threads), "--runMode", "genomeGenerate", "--genomeDir", params.ref_dir, "--genomeFastaFiles", input.ref, "--genomeSAindexNbases", "6"]
            else:
                raise NotImplementedError(f"Mapper '{mapper}' is not supported for indexing.")
            subprocess.run(cmd, check=True)
            logger.info(f"Index for {mapper} created successfully.")
        except Exception as e:
            logger.error(f"Failed to create index for {mapper}: {e}")
            raise e

rule map_reads:
    """
    Maps trimmed reads to the reference. For bowtie2, it runs both --local and
    --end-to-end modes and merges the results.
    """
    input:
        r1="results/03-trim_data/{species}/{sra}_R1.trimmed.fastq.gz",
        r2="results/03-trim_data/{species}/{sra}_R2.trimmed.fastq.gz",
        idx_done="results/06-mapping/index/{reference}_{mapper}/index.done"
    output:
        bam="results/06-mapping/bams/{species}/{sra}/{reference}/{mapper}/mapped.bam"
    threads: 16
    run:
        log_file = f"results/06-mapping/logs/{wildcards.species}_{wildcards.sra}_{wildcards.reference}_{wildcards.mapper}.log"
        logger = setup_logger(f"map_reads_{wildcards.sra}_{wildcards.reference}_{wildcards.mapper}", log_file)
        try:
            mapper = wildcards.mapper
            reference = wildcards.reference
            ref_prefix = f"results/06-mapping/index/{reference}_{mapper}/{reference}"
            star_index_dir = f"results/06-mapping/index/{reference}_{mapper}/"
            
            logger.info(f"Mapping reads for {wildcards.sra} to {reference} using {mapper}.")
            map_cmd = ""
            
            with tempfile.TemporaryDirectory() as tmp_dir:
                if mapper == "bowtie2":
                    logger.info("Bowtie2 mode: running --local and --end-to-end mappings.")
                    local_bam = os.path.join(tmp_dir, "local.bam")
                    e2e_bam = os.path.join(tmp_dir, "end_to_end.bam")

                    local_cmd = f"bowtie2 --local -p {threads} -x {ref_prefix} -1 {input.r1} -2 {input.r2} | samtools view -bS - > {local_bam}"
                    e2e_cmd = f"bowtie2 --end-to-end -p {threads} -x {ref_prefix} -1 {input.r1} -2 {input.r2} | samtools view -bS - > {e2e_bam}"
                    subprocess.run(local_cmd, shell=True, check=True, executable='/bin/bash')
                    subprocess.run(e2e_cmd, shell=True, check=True, executable='/bin/bash')
                    
                    logger.info("Merging bowtie2 results.")
                    map_cmd = f"samtools merge -u -@ {threads} - {local_bam} {e2e_bam}"

                elif mapper == "bwa":
                    map_cmd = f"bwa mem -t {threads} {ref_prefix} {input.r1} {input.r2}"
                
                elif mapper == "star":
                    star_prefix = os.path.join(tmp_dir, "star_")
                    star_run_cmd = f"STAR --runThreadN {threads} --genomeDir {star_index_dir} --readFilesIn {input.r1} {input.r2} --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix {star_prefix}"
                    subprocess.run(star_run_cmd, shell=True, check=True)
                    map_cmd = f"cat {star_prefix}Aligned.out.bam"

                logger.info("Filtering for mapped reads, sorting, and writing final BAM.")
                full_command = (
                    f"{map_cmd} | "
                    f"samtools view -@ {threads} -F 4 -b - | "
                    f"samtools sort -@ {threads} -o {output.bam}"
                )
                subprocess.run(full_command, shell=True, check=True, executable='/bin/bash')
                subprocess.run(["samtools", "index", output.bam], check=True)

            logger.info(f"Successfully created final mapped BAM for {wildcards.sra}.")
        except Exception as e:
            logger.error(f"Mapping failed for {wildcards.sra}: {e}")
            raise e

rule assemble_contigs:
    """
    Extracts mapped reads from the BAM file and uses them for de novo assembly.
    This version now uses a simplified approach, focusing only on paired-end reads
    and ensuring the Trinity output path is valid inside the container.
    """
    input:
        bam="results/06-mapping/bams/{species}/{sra}/{reference}/{mapper}/mapped.bam",
        container_sif=lambda wildcards: os.path.join(
            config["container_dir"], f"{get_mapper_and_reference(wildcards)['assembler']}.sif"
        ) if get_mapper_and_reference(wildcards).get('assembler') in config["containers"] else []
    output:
        contigs="results/07-assembly/{species}/{sra}/{reference}/{mapper}/{assembler}/contigs.fasta"
    params:
        outdir="results/07-assembly/{species}/{sra}/{reference}/{mapper}/{assembler}",
        assembler="{assembler}",
        mitoz_clade=config.get("params", {}).get("mitoz", {}).get("clade", "Chordata"),
        mitoz_gcode=config.get("params", {}).get("mitoz", {}).get("genetic_code", 2)
    threads: 16
    shell:
        """
        set -e
        
        TMP_DIR=$(mktemp -d)
        R1=$TMP_DIR/R1.fastq.gz
        R2=$TMP_DIR/R2.fastq.gz
        # Removido o suporte a singletons por enquanto para simplificar

        echo "Extracting mapped PAIRED-END reads from {input.bam}..."
        # -f 1: Extrai apenas reads que estão devidamente pareadas.
        # Removido o -s para ignorar singletons.
        samtools fastq -F 4 -f 1 -@ {threads} -1 $R1 -2 $R2 {input.bam}

        if [ ! -s "$R1" ]; then
            echo "No mapped paired-end reads found. Creating an empty contig file."
            echo ">no_reads_found_for_{wildcards.sra}" > {output.contigs}
            rm -r $TMP_DIR
            exit 0
        fi

        echo "Reads extracted. Running {wildcards.assembler} assembler..."
        
        # O diretório de saída no host é montado em um caminho compatível com o Trinity dentro do contêiner
        SINGULARITY_CMD="sudo singularity exec -e -B $TMP_DIR:/data/input -B {params.outdir}:/data/trinity_output {input.container_sif}"

        if [ "{wildcards.assembler}" = "trinity" ]; then
            # O --output agora aponta para um caminho válido dentro do contêiner
            $SINGULARITY_CMD Trinity --seqType fq --left /data/input/R1.fastq.gz --right /data/input/R2.fastq.gz --max_memory 50G --CPU {threads} --output /data/trinity_output
            mv {params.outdir}/Trinity.fasta {output.contigs}
        
        elif [ "{wildcards.assembler}" = "mitoz" ]; then
            $SINGULARITY_CMD MitoZ.py assemble --genetic_code {params.mitoz_gcode} --clade {params.mitoz_clade} --thread_number {threads} --outprefix /data/trinity_output/mitoz_assembly --fastq1 /data/input/R1.fastq.gz --fastq2 /data/input/R2.fastq.gz
            mv {params.outdir}/mitoz_assembly.result/work71.mitogenome.fa {output.contigs}
        
        else # Assemblers que não precisam de contêiner
            spades.py --careful -1 $R1 -2 $R2 -o {params.outdir} -t {threads}
            mv {params.outdir}/contigs.fasta {output.contigs}
        fi

        rm -r $TMP_DIR
        echo "Assembly for {wildcards.sra} completed successfully."
        """