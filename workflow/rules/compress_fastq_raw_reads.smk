import os
import shutil
import subprocess

rule compress_fastq_raw_reads:
    """Compresses the downloaded FASTQ files using gzip or pigz."""
    input:
        r1="results/01-raw_data/{species}/{sra}_R1.fastq",
        r2="results/01-raw_data/{species}/{sra}_R2.fastq"
    output:
        r1_gz="results/01-raw_data/{species}/{sra}_R1.fastq.gz",
        r2_gz="results/01-raw_data/{species}/{sra}_R2.fastq.gz"
    threads: 4
    #benchmark:
    #    repeat("benchmarks/{species}/{sra}.compress_fastq_raw_reads.tsv", 5)
    run:
        logger = setup_logger(f"compress_fastq_raw_reads_{wildcards.sra}")

        #if not os.path.exists("benchmarks"):
        #    logger.info("Creating benchmarks directory.")
        #    os.makedirs("benchmarks", exist_ok=True, mode=0o755)
        
        try:
            logger.info(f"Starting compression for SRA accession '{wildcards.sra}'.")
            
            if os.path.exists(str(output.r1_gz)) and os.path.exists(str(output.r2_gz)):
                logger.info("Compressed FASTQ files already exist. Skipping compression.")
            else:
                logger.info(f"Compressing FASTQ files using {threads} threads...")
                
                if shutil.which("pigz"):
                    
                    with open(str(output.r1_gz), "wb") as f1:
                        run = subprocess.run(["pigz", "-9", "-p", str(threads), "-c", str(input.r1)], 
                                     check=True, stdout=f1, stderr=subprocess.PIPE)
                        if run.returncode != 0:
                            raise RuntimeError(f"pigz failed with error: {run.stderr.decode('utf-8')}")
                    
                    with open(str(output.r2_gz), "wb") as f2:
                        run = subprocess.run(["pigz", "-9", "-p", str(threads), "-c", str(input.r2)], 
                                     check=True, stdout=f2, stderr=subprocess.PIPE)
                        if run.returncode != 0:
                            raise RuntimeError(f"pigz failed with error: {run.stderr.decode('utf-8')}")
                
                elif shutil.which("gzip"):
                    
                    with open(str(output.r1_gz), "wb") as f1:
                        run = subprocess.run(["gzip", "-9", "-c", str(input.r1)], 
                                     check=True, stdout=f1, stderr=subprocess.PIPE)
                        if run.returncode != 0:
                            raise RuntimeError(f"gzip failed with error: {run.stderr.decode('utf-8')}")

                    with open(str(output.r2_gz), "wb") as f2:
                        run = subprocess.run(["gzip", "-9", "-c", str(input.r2)], 
                                     check=True, stdout=f2, stderr=subprocess.PIPE)
                        if run.returncode != 0:
                            raise RuntimeError(f"gzip failed with error: {run.stderr.decode('utf-8')}")
                else:
                    raise EnvironmentError("Neither pigz nor gzip is available for compression.")

            
            logger.info(f"Compression for SRA accession '{wildcards.sra}' finished successfully.")
        except Exception as e:
            logger.error(f"Compression job for SRA accession '{wildcards.sra}' failed: {e}")
            
            if os.path.exists(str(output.r1_gz)):
                shutil.remove(str(output.r1_gz))
            if os.path.exists(str(output.r2_gz)):
                shutil.remove(str(output.r2_gz))

            raise e