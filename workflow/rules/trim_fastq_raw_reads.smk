rule trim_fastq_raw_reads:
    """Uses fastp to perform adapter trimming and quality filtering."""
    input:
        r1="results/01-raw_data/{species}/{sra}_R1.fastq.gz",
        r2="results/01-raw_data/{species}/{sra}_R2.fastq.gz"
    output:
        r1_trim="results/03-trim_data/{species}/{sra}_R1.trimmed.fastq.gz",
        r2_trim="results/03-trim_data/{species}/{sra}_R2.trimmed.fastq.gz",
        json="results/03-trim_data/{species}/fastp/{sra}.fastp.json",
        html="results/03-trim_data/{species}/fastp/{sra}.fastp.html"
    params:
        extra=config.get("params", {}).get("fastp", "")
    threads: 8
    #benchmark:
    #    repeat("benchmarks/{species}/{sra}.trim_fastq_raw_reads.tsv", 5)
    run:
        logger = setup_logger(f"trim_fastq_raw_reads_{wildcards.sra}")

        #if not os.path.exists("benchmarks"):
        #    logger.info("Creating benchmarks directory.")
        #    os.makedirs("benchmarks", exist_ok=True, mode=0o755)

        try:
            logger.info(f"Running fastp for {wildcards.sra} with {threads} threads.")
            
            if not os.path.exists(os.path.dirname(str(output.r1_trim))):
                os.makedirs(os.path.dirname(str(output.r1_trim)), exist_ok=True, mode=0o755)
            
            cmd = [
                "fastp",
                "--in1", input.r1,
                "--in2", input.r2,
                "--out1", output.r1_trim,
                "--out2", output.r2_trim,
                "--json", output.json,
                "--html", output.html,
                "--thread", str(threads),
                "--detect_adapter_for_pe",
                "--trim_poly_g",
                "--trim_poly_x",
                ]
            
            if isinstance(params.extra, dict):
                for k, v in params.extra.items():
                    cmd.append(str(k))
                    if v is not True:
                        cmd.append(str(v))
            elif isinstance(params.extra, str) and params.extra:
                cmd.extend(params.extra.split())
            
            run = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if run.returncode != 0:
                raise RuntimeError(f"fastp failed with error: {run.stderr.decode('utf-8')}")

            time.sleep(5)
            logger.info(f"fastp for {wildcards.sra} completed successfully.")
        
        except subprocess.CalledProcessError as e:
            logger.error(f"fastp for {wildcards.sra} failed:\n{e.stderr}")
            raise e