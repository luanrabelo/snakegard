import os
import shutil
import subprocess
import time

rule index_reference:
    """Creates an index of the reference sequence for a specific mapper."""
    input:
        ref="results/00-references/{reference}.fasta"
    output:
        done="results/06-mapping/{species}/{sra}/index/{reference}_{mapper}/index.done"
    params:
        ref_prefix=lambda wildcards: f"results/06-mapping/{wildcards.species}/{wildcards.sra}/index/{wildcards.reference}_{wildcards.mapper}/{wildcards.reference}",
        ref_dir=lambda wildcards: f"results/06-mapping/{wildcards.species}/{wildcards.sra}/index/{wildcards.reference}_{wildcards.mapper}/"
    threads: 8
    run:
        logger = setup_logger(f"index_{wildcards.mapper}_{wildcards.reference}")
        try:
            mapper = wildcards.mapper
            
            logger.info(f"Creating index for {input.ref} using {mapper}.")
            
            # Criar diretório usando o parâmetro correto
            ref_dir = params.ref_dir if callable(params.ref_dir) else params.ref_dir
            if callable(ref_dir):
                ref_dir = ref_dir(wildcards)
            
            if not os.path.exists(ref_dir):
                os.makedirs(ref_dir, exist_ok=True, mode=0o755)
            
            # Verificar se o mapper está disponível
            if mapper == "bowtie2":
                if not shutil.which("bowtie2-build"):
                    raise EnvironmentError("bowtie2-build is not available. Please install bowtie2: conda install -c bioconda bowtie2")
                
                ref_prefix = params.ref_prefix if callable(params.ref_prefix) else params.ref_prefix
                if callable(ref_prefix):
                    ref_prefix = ref_prefix(wildcards)
                
                cmd = [
                    "bowtie2-build",
                    "--threads",
                    str(threads),
                    str(input.ref),
                    str(ref_prefix)
                    ]
                
            elif mapper == "bwa":
                if not shutil.which("bwa"):
                    raise EnvironmentError("bwa is not available. Please install bwa: conda install -c bioconda bwa")
                
                ref_prefix = params.ref_prefix if callable(params.ref_prefix) else params.ref_prefix
                if callable(ref_prefix):
                    ref_prefix = ref_prefix(wildcards)
                
                cmd = [
                    "bwa",
                    "index",
                    "-p",
                    str(ref_prefix),
                    str(input.ref)
                    ]
                
            elif mapper == "star":
                if not shutil.which("STAR"):
                    raise EnvironmentError("STAR is not available. Please install star: conda install -c bioconda star")
                
                cmd = [
                    "STAR",
                    "--runThreadN",
                    str(threads),
                    "--runMode",
                    "genomeGenerate",
                    "--genomeDir",
                    str(ref_dir),
                    "--genomeFastaFiles",
                    str(input.ref),
                    "--genomeSAindexNbases",
                    "6"
                    ]
            else:
                raise NotImplementedError(f"Mapper '{mapper}' is not supported for indexing.")
            
            logger.info(f"Running with tool {mapper} and command: {' '.join(cmd)}")
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            if result.returncode != 0:
                logger.error(f"Indexing failed with error: {result.stderr}")
                raise RuntimeError(f"Indexing failed for {mapper} on reference {wildcards.reference}")
            
            # Criar arquivo de confirmação
            with open(output.done, "w") as f:
                f.write(f"Index created for {wildcards.reference} using {mapper} on {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"Command: {' '.join(cmd)}\n")
                f.write(f"Output directory: {ref_dir}\n")
            
            logger.info(f"Index for {mapper} created successfully.")
        
        except Exception as e:
            logger.error(f"Failed to create index for {mapper}: {e}")
            if os.path.exists(output.done):
                os.remove(output.done)
            raise e