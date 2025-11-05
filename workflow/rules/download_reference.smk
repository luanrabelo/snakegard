import os
import shutil
import subprocess

rule download_reference:
    """Downloads a nucleotide sequence from NCBI using its accession number."""
    output:
        fasta="results/00-references/{reference}.fasta"
    params:
        db="nucleotide"
    threads: 2
    #benchmark:
    #    repeat("benchmarks/{species}/{sra}.download_reference.tsv", 5)
    run:
        logger = setup_logger(f"download_reference_{wildcards.reference}")
        
        try:
            logger.info(f"Starting job for reference {wildcards.reference}.")

             #if not os.path.exists("benchmarks"):
             #  logger.info("Creating benchmarks directory.")
             #  os.makedirs("benchmarks", exist_ok=True, mode=0o755)
            
            if os.path.exists(str(output.fasta)):
                logger.info("Reference file already exists. Skipping.")
            else:
                max_retries = 5
                retry_wait_seconds = 20
                downloaded_content = ""
                
                for attempt in range(max_retries):
                    logger.info(f"Downloading Reference Sequence {wildcards.reference} (Attempt {attempt + 1}/{max_retries})...")
                    
                    try:
                        cmd = [
                            "efetch",
                            "-db",
                            str(params.db),
                            "-id",
                            str(wildcards.reference),
                            "-format",
                            "fasta"
                            ]
                        
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