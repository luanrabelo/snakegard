<p align="center">
  <img src="assets/SnakeGARD.png" alt="SnakeGARD Logo" width="100%">
</p>

<p align="center">
  <a href="https://www.buymeacoffee.com/lprabelo" target="_blank">
    <img src="https://img.buymeacoffee.com/button-api/?text=Buy me a coffee&emoji=☕&slug=lprabelo&button_colour=f1f1f1&font_colour=000000&font_family=Lato&outline_colour=000000&coffee_colour=000000" />
  </a>
</p>

# Contents Overview
- [System Overview](#system-overview)
  - [Pipeline Workflow](#pipeline-workflow)
  - [Features](#features)
- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Setting up the Conda Environment](#setting-up-the-conda-environment)
  - [Singularity/Apptainer Setup](#singularityapptainer-setup)
- [Usage](#usage)
  - [Configuration](#configuration)
  - [Running the Pipeline](#running-the-pipeline)
  - [Resuming Failed Runs](#resuming-failed-runs)
  - [Pipeline Output](#pipeline-output)
- [Customization](#customization)
  - [Advanced Configuration Options](#advanced-configuration-options)
  - [Using Different Assemblers](#using-different-assemblers)
- [Troubleshooting](#troubleshooting)
  - [Common Issues](#common-issues)
- [Licence](#licence)
- [Contact](#contact)

# System Overview
##### [:rocket: Go to Contents Overview](#contents-overview)

**SnakeGARD** (Snake Genome Assembly from RNA-Seq Data) is a Snakemake-based workflow designed for automated and scalable assembly of mitochondrial genomes from RNA-Seq data. The pipeline integrates various bioinformatics tools to provide a robust and efficient solution for mitochondrial genome assembly.

## Pipeline Workflow

SnakeGARD follows a comprehensive workflow that includes:

1. **Data Acquisition**: Automatic download of reference sequences and RNA-Seq data from NCBI
2. **Quality Control**: FastQC analysis of raw and trimmed reads
3. **Read Preprocessing**: Adapter trimming and quality filtering using fastp
4. **Mapping**: Alignment of trimmed reads to reference sequences using Bowtie2, BWA, or STAR
5. **Assembly**: De novo assembly of mapped reads using various assemblers (SPAdes, Trinity, MitoZ)
6. **Quality Assessment**: Comprehensive quality reporting with MultiQC

## Features

- **Fully Automated**: From data download to assembly with minimal user intervention
- **Containerized**: Uses Singularity/Apptainer containers for reproducibility
- **Flexible**: Supports multiple mappers and assemblers
- **Scalable**: Efficiently processes multiple samples in parallel
- **Resumable**: Can continue from the last successful step after failures

# Installation
##### [:rocket: Go to Contents Overview](#contents-overview)

## Prerequisites

Before installing SnakeGARD, ensure you have the following software installed:

- **Git**: For cloning the repository
- **Conda/Miniconda**: For managing the environment and dependencies
- **Singularity/Apptainer**: For running containerized applications

## Setting up the Conda Environment

1. **Clone the repository**:
   ```bash
   git clone https://github.com/luanrabelo/snakegard.git
   cd snakegard
   ```

2. **Create the Conda environment**:
   ```bash
   conda env create -f environment.yml
   ```

3. **Activate the environment**:
   ```bash
   conda activate snakegard
   ```

## Singularity/Apptainer Setup

SnakeGARD uses Singularity/Apptainer containers to ensure reproducibility. The containers are automatically downloaded during the workflow execution, but you need to have Singularity/Apptainer installed on your system.

For Ubuntu/Debian:
```bash
sudo apt update
sudo apt install -y singularity-container
```

For other operating systems, please refer to the [Singularity installation guide](https://sylabs.io/guides/latest/user-guide/quick_start.html).

# Usage
##### [:rocket: Go to Contents Overview](#contents-overview)

## Configuration

SnakeGARD is configured through the `config.yaml` file. This file contains parameters for data sources, reference sequences, mappers, and assemblers.

### Sample Configuration

```yaml
container_dir: "workflow/containers"

containers:
  trinity: "https://data.broadinstitute.org/Trinity/TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg"
  mitoz: "https://www.dropbox.com/scl/fo/40jgsbnoocqr9qyns7tql/AB3QgelqF1Uex_90OSwR9MA?rlkey=mw03ho8uwwomr7amw1ru430o1&e=1&dl=1"

params:
  fastp:
    length_required: 100
    qualified_quality_phred: 25
  mitoz:
    clade: "Chordata"
    genetic_code: 2

samples:
  - species: "Anaxyrus_microscaphus"
    sra:
      - id: "SRR18432214"
        reference_type: "mitochondrion"
        reference: "NC_047224.1"
        mapper: "bowtie2"
        assembler: "mitoz"
```

### Configuration Parameters

- **container_dir**: Directory where container files (.sif) will be stored
- **containers**: URLs for downloading assembler containers
- **params**: Tool-specific parameters
  - **fastp**: Read trimming parameters
  - **mitoz**: MitoZ assembler parameters
- **samples**: List of samples to process
  - **species**: Species name (use underscores instead of spaces)
  - **sra**: List of SRA entries
    - **id**: SRA accession number
    - **reference_type**: Type of reference sequence
    - **reference**: Reference sequence accession number
    - **mapper**: Mapping tool to use (bowtie2, bwa, or star)
    - **assembler**: Assembly tool to use (spades, rnaspades, mitoz, or trinity)

## Running the Pipeline

To run SnakeGARD with default parameters:

```bash
snakemake --use-conda --use-singularity --cores <number_of_cores>
```

For a dry run to check the workflow without executing any commands:

```bash
snakemake -n
```

To generate a workflow graph:

```bash
snakemake --dag | dot -Tpng > dag.png
```

## Resuming Failed Runs

If the pipeline fails, you can resume from the last successful step:

```bash
snakemake --use-conda --use-singularity --cores <number_of_cores> --rerun-incomplete
```

## Pipeline Output

SnakeGARD organizes results in the `results/` directory with the following structure:

- **00-references/**: Reference sequences downloaded from NCBI
- **01-raw_data/**: Raw FASTQ files from SRA
- **02-fastqc_raw/**: FastQC reports for raw reads
- **03-trim_data/**: Trimmed FASTQ files
- **04-fastqc_trimmed/**: FastQC reports for trimmed reads
- **05-multiqc/**: Aggregated quality control reports
- **06-mapping/**: Mapping results (BAM files)
- **07-assembly/**: Final assembly results (FASTA files)

# Customization
##### [:rocket: Go to Contents Overview](#contents-overview)

## Advanced Configuration Options

You can customize SnakeGARD by modifying the following:

### Resource Allocation

Control the number of threads used by each rule in the Snakefile:

```yaml
# Add to config.yaml
resources:
  mapping_threads: 16
  assembly_threads: 32
```

### Custom Containers

Add custom containers for additional assemblers:

```yaml
# Add to config.yaml
containers:
  custom_assembler: "https://path/to/your/container.sif"
```

## Using Different Assemblers

SnakeGARD supports multiple assemblers that can be specified in the `config.yaml`:

- **SPAdes**: General-purpose assembler for small genomes
- **RNASPAdes**: Assembler optimized for transcriptome data
- **MitoZ**: Specialized for mitochondrial genome assembly
- **Trinity**: Transcriptome assembler

# Troubleshooting
##### [:rocket: Go to Contents Overview](#contents-overview)

## Common Issues

### Missing SRA Tools

If you encounter errors with SRA data download:

```
Error: fasterq-dump command not found
```

Solution:
```bash
conda install -c bioconda sra-tools
```

### Singularity Permission Issues

If you encounter permission errors with Singularity:

```
ERROR: Failed to create user namespace
```

Solution:
```bash
# For Linux systems
sudo singularity config --set allow_root_capabilities=true
# Or run with sudo
sudo snakemake --use-conda --use-singularity --cores <number_of_cores>
```

### Out of Memory Errors

If assemblies fail due to memory issues:

```
Error: Process killed (out of memory)
```

Solution: Adjust the memory limit in your Snakemake command:

```bash
snakemake --use-conda --use-singularity --cores <number_of_cores> --resources mem_mb=<memory_in_MB>
```

## Licence
##### [:rocket: Go to Contents Overview](#contents-overview)
***SnakeGARD*** is released under the **MIT License**. This license permits reuse within proprietary software provided that all copies of the licensed software include a copy of the MIT License terms and the copyright notice.
***

# Contact
##### [:rocket: Go to Contents Overview](#contents-overview)
For reporting bugs, requesting assistance, or providing feedback, please reach out to **Luan Rabelo**:
```
luanrabelo@outlook.com
```
***