# assembly: rules for long read assembly and end repair.
# Copyright Jackson M. Tsuji and Lee H. Bergstrand, 2023

import os

# SAMPLE_NAMES is instantiated in rotary.smk

rule assembly_flye:
    input:
        "{sample}/qc/{sample}_qc_long.fastq.gz"
    output:
        assembly="{sample}/assembly/flye/{sample}_flye_assembly.fasta",
        info="{sample}/assembly/flye/{sample}_assembly_info.txt",
        output_dir=directory("{sample}/assembly/flye"),
        assembly_dir=temp(directory("{sample}/assembly/flye/00-assembly")),
        consensus_dir=temp(directory("{sample}/assembly/flye/10-consensus")),
        repeat_dir=temp(directory("{sample}/assembly/flye/20-repeat")),
        contigger_dir=temp(directory("{sample}/assembly/flye/30-contigger")),
        polishing_dir=temp(directory("{sample}/assembly/flye/40-polishing")) if config.get("flye_polishing_rounds") > 0 else []
    conda:
        "../envs/assembly_flye.yaml"
    log:
        "{sample}/logs/assembly/assembly_flye.log"
    benchmark:
        "{sample}/benchmarks/assembly/assembly_flye.benchmark.txt"
    params:
        input_mode=config.get("flye_input_mode"),
        read_error="" if config.get("flye_read_error") == "auto" else "--read-error " + config.get("flye_read_error"),
        meta_mode="--meta" if config.get("flye_meta_mode") == "True" else "",
        polishing_rounds=config.get("flye_polishing_rounds")
    threads:
        config.get("threads",1)
    shell:
        """
        flye --{params.input_mode} {input} {params.read_error} --out-dir {output.output_dir} {params.meta_mode} \
          --iterations {params.polishing_rounds} -t {threads} > {log} 2>&1
        mv {output.output_dir}/assembly.fasta {output.assembly} 
        mv {output.output_dir}/assembly_info.txt {output.info} 
        """


# TODO - eventually make this rule more generic so that outputs from other assemblers can go to the same output files
# TODO - the math for memory per thread should really be done somewhere other than 'resources' - what is best practice?
rule assembly_end_repair:
    input:
        qc_long_reads="{sample}/qc/{sample}_qc_long.fastq.gz",
        assembly="{sample}/assembly/flye/{sample}_flye_assembly.fasta",
        info="{sample}/assembly/flye/{sample}_assembly_info.txt"
    output:
        assembly="{sample}/assembly/end_repair/{sample}_end_repair_assembly.fasta",
        info="{sample}/assembly/end_repair/{sample}_repaired_info.tsv",
        output_dir=directory("{sample}/assembly/end_repair")
    conda:
        "../envs/spokewrench.yaml"
    log:
        "{sample}/logs/assembly/end_repair.log"
    benchmark:
        "{sample}/benchmarks/assembly/end_repair.txt"
    params:
        flye_input_mode=config.get("flye_input_mode"),
        flye_read_error="0" if config.get("flye_read_error") == "auto" else config.get("flye_read_error"),
        min_id=config.get("circlator_merge_min_id"),
        min_length=config.get("circlator_merge_min_length"),
        ref_end=config.get("circlator_merge_ref_end"),
        reassemble_end=config.get("circlator_merge_reassemble_end"),
        keep_going="--keep_going_with_failed_contigs" if config.get("keep_unrepaired_contigs") == "True" else "",
        mem_per_thread=int(config.get("memory") /config.get("threads",1))
    threads:
        config.get("threads",1)
    resources:
        mem=config.get("memory")
    shell:
        """
        spokewrench repair --long_read_filepath {input.qc_long_reads} \
          --assembly_fasta_filepath {input.assembly} \
          --assembly_info_filepath {input.info} \
          --output_dir {output.output_dir} \
          --flye_read_mode {params.flye_input_mode} \
          --flye_read_error {params.flye_read_error} \
          --circlator_min_id {params.min_id} \
          --circlator_min_length {params.min_length} \
          --circlator_ref_end {params.ref_end} \
          --circlator_reassemble_end {params.reassemble_end} \
          --threads {threads} \
          --threads_mem {params.mem_per_thread} \
          --overwrite \
          {params.keep_going} \
          > {log} 2>&1
        mv {output.output_dir}/repaired.fasta {output.assembly} 
        mv {output.output_dir}/repaired_info.tsv {output.info}
        """

rule run_quast:
    input:
        flye_assembly = "{sample}/assembly/flye/{sample}_flye_assembly.fasta"
    output:
        stats_dir=directory("{sample}/assembly/stats/"),
        report=multiext("{sample}/assembly/stats/report", '.html', '.tex', '.tsv', '.txt', '.pdf'),
        transposed_report=multiext("{sample}/assembly/stats/transposed_report",'.tex', '.tsv', '.txt')
    conda:
        "../envs/qc.yaml"
    log:
        "{sample}/logs/assembly/assembly_stats.log"
    benchmark:
        "{sample}/benchmarks/assembly/assembly_stats.benchmark.txt"
    threads:
        2
    shell:
        """
        quast -t {threads} --space-efficient \
        --labels "{wildcards.sample}_flye, {wildcards.sample}_end_repair" \
        -o {output.stats_dir} {input} > {log} 2>&1
        """


rule finalize_assembly:
    input:
        assembly="{sample}/assembly/end_repair/{sample}_end_repair_assembly.fasta",
        info="{sample}/assembly/end_repair/{sample}_repaired_info.tsv",
        stats="{sample}/assembly/stats/report.pdf"
    output:
        assembly="{sample}/assembly/{sample}_assembly.fasta",
        info="{sample}/assembly/{sample}_circular_info.tsv"
    run:
        source_relpath = os.path.relpath(str(input.assembly),os.path.dirname(str(output.assembly)))
        os.symlink(source_relpath,str(output.assembly))

        source_relpath = os.path.relpath(str(input.info),os.path.dirname(str(output.info)))
        os.symlink(source_relpath,str(output.info))


rule aggregate_assembly_stats_multiqc:
    input:
        expand("{sample}/assembly/stats/",sample=SAMPLE_NAMES)
    output:
        multqc_report="stats/assembly/all_samples_assembly_stats_multiqc_report.html",
        multqc_report_data="stats/assembly/all_samples_assembly_stats_multiqc_report_data.zip"
    conda:
        "../envs/qc.yaml"
    log:
        "logs/assembly/assembly_stats_multiqc.log"
    benchmark:
        "benchmarks/assembly/assembly_stats_multiqc.txt"
    params:
        qc_stats_dir="stats/assembly/"
    shell:
        """
        multiqc --outdir {params.qc_stats_dir} \
        --title 'all_samples_assembly_stats' \
        --data-format tsv \
        --zip-data-dir {input} > {log} 2>&1
        """


rule assembly:
    input:
        assembly_fasta = expand("{sample}/assembly/{sample}_assembly.fasta",sample=SAMPLE_NAMES),
        assembly_circular_info = expand("{sample}/assembly/{sample}_circular_info.tsv", sample=SAMPLE_NAMES),
        assembly_multqc_report="stats/assembly/all_samples_assembly_stats_multiqc_report.html",
        assembly_multqc_report_data="stats/assembly/all_samples_assembly_stats_multiqc_report_data.zip"
    output:
        temp(touch("checkpoints/assembly"))
