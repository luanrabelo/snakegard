import os
import subprocess

rule fastqc_trimmed_reads:
    """Runs FastQC on the trimmed, compressed FASTQ files."""
    input:
        r1="results/03-trim_data/{species}/{sra}_R1.trimmed.fastq.gz",
        r2="results/03-trim_data/{species}/{sra}_R2.trimmed.fastq.gz"
    output:
        r1_zip="results/04-fastqc_trimmed/{species}/{sra}_R1.trimmed_fastqc.zip",
        r2_zip="results/04-fastqc_trimmed/{species}/{sra}_R2.trimmed_fastqc.zip",
        r1_html="results/04-fastqc_trimmed/{species}/{sra}_R1.trimmed_fastqc.html",
        r2_html="results/04-fastqc_trimmed/{species}/{sra}_R2.trimmed_fastqc.html"
    params:
        outdir="results/04-fastqc_trimmed/{species}"
    threads: 4
    #benchmark:
    #    repeat("benchmarks/{species}/{sra}.fastqc_trimmed_reads.tsv", 5)
    run:
        logger = setup_logger(f"fastqc_trimmed_{wildcards.sra}")
        
        #if not os.path.exists("benchmarks"):
        #    logger.info("Creating benchmarks directory.")
        #    os.makedirs("benchmarks", exist_ok=True, mode=0o755)

        try:
            logger.info(f"Running FastQC on trimmed reads for {wildcards.sra}.")

            if not os.path.exists(params.outdir):
                os.makedirs(name=params.outdir, exist_ok=True, mode=0o755)
            
            cmd = [
                "fastqc",
                "--threads",
                str(threads),
                "--outdir",
                str(params.outdir),
                str(input.r1),
                str(input.r2)
                ]
            
            run = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            if run.returncode != 0:
                raise RuntimeError(f"FastQC failed with error: {run.stderr.decode('utf-8')}")
            
            time.sleep(10)
            logger.info(f"FastQC on trimmed reads for {wildcards.sra} completed successfully.")

        except subprocess.CalledProcessError as e:
            logger.error(f"FastQC on trimmed reads for {wildcards.sra} failed: {e.stderr.decode()}")
            raise e