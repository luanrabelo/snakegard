import os
import shutil
import subprocess
import tempfile

rule download_fastq_raw_reads:
    """Downloads paired-end FASTQ data from SRA for a given accession number."""
    output:
        r1=temp("results/01-raw_data/{species}/{sra}_R1.fastq"),
        r2=temp("results/01-raw_data/{species}/{sra}_R2.fastq")
    threads: 8
    #benchmark:
    #    repeat("benchmarks/{species}/{sra}.download_fastq_raw_reads.tsv", 5)
    run:
        logger = setup_logger(f"download_fastq_raw_reads_{wildcards.sra}")

        #if not os.path.exists("benchmarks"):
        #    logger.info("Creating benchmarks directory.")
        #    os.makedirs("benchmarks", exist_ok=True, mode=0o755)

        try:
            logger.info(f"Starting job for SRA accession {wildcards.sra}.")
            
            if os.path.exists(str(output.r1)) and os.path.exists(str(output.r2)):
                 logger.info("Temporary FASTQ files already exist. Skipping download.")
            else:
                os.makedirs(f"results/01-raw_data/{wildcards.species}", exist_ok=True, mode=0o755)
                
                if shutil.which("fasterq-dump") is None:
                    raise EnvironmentError("fasterq-dump is not available in the system PATH.")
                else:
                    with tempfile.TemporaryDirectory(suffix=".tmp", prefix=f"snakeseek_", dir=f"{os.path.join(os.getcwd(), 'tmp')}") as tmp_dir:
                        logger.info(f"Downloading {wildcards.sra} to temporary directory {tmp_dir} using {threads} threads...")
                        
                        cmd = [
                            "fasterq-dump",
                            str(wildcards.sra),
                            "--split-files",
                            "--threads",
                            str(threads),
                            "--outdir",
                            str(tmp_dir)
                            ]
                        
                        run = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=tmp_dir)

                        if run.returncode != 0:
                            raise RuntimeError(f"fasterq-dump failed with error: {run.stderr.decode('utf-8')}")
                        
                        src_r1 = os.path.join(tmp_dir, f"{wildcards.sra}_1.fastq")
                        src_r2 = os.path.join(tmp_dir, f"{wildcards.sra}_2.fastq")
                        
                        if not os.path.exists(src_r1) or not os.path.exists(src_r2):
                            raise FileNotFoundError(f"fasterq-dump did not produce the expected FASTQ files in {tmp_dir}")
                        
                        logger.info("Download complete. Moving files to final destination.")
                        
                        shutil.move(src_r1, str(output.r1))
                        shutil.move(src_r2, str(output.r2))
                
                logger.info(f"Job for SRA accession {wildcards.sra} finished successfully.")
        
        except Exception as e:
            logger.error(f"Job for SRA accession {wildcards.sra} failed: {e}")
            
            if os.path.exists(str(output.r1)):
                shutil.remove(str(output.r1))
            if os.path.exists(str(output.r2)):
                shutil.remove(str(output.r2))
            
            raise e