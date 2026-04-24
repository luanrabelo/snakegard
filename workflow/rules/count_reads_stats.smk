rule count_reads_stats:
    input:
        fastp_json="results/03-trim_data/{species}/fastp/{sra}.fastp.json",
        bam="results/06-mapping/{species}/{sra}/bam/{species}_{reference}_{mapper}_mapped.bam"
    output:
        tsv="results/08-stats/{species}_{sra}_{reference}_{mapper}_counts.tsv"
    run:
        import json
        import subprocess

        with open(input.fastp_json) as f:
            fastp_data = json.load(f)
            raw_reads = fastp_data['summary']['before_filtering']['total_reads']
            trimmed_reads = fastp_data['summary']['after_filtering']['total_reads']

        cmd = f"samtools view -c -F 4 {input.bam}"
        mapped_reads = subprocess.check_output(cmd, shell=True).decode().strip()

        with open(output.tsv, 'w') as f:
            f.write("Species\tSRA\tReference\tMapper\tRaw_Reads\tTrimmed_Reads\tMapped_Reads\n")
            f.write(f"{wildcards.species}\t{wildcards.sra}\t{wildcards.reference}\t{wildcards.mapper}\t{raw_reads}\t{trimmed_reads}\t{mapped_reads}\n")