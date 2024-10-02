# annotation: rules for genome annotation.
# Copyright Jackson M. Tsuji and Lee H. Bergstrand, 2023

import os
from pungi.utils import is_config_parameter_true
from rotary.annotation import AnnotationMap, combine_tabular_reports

VERSION_DFAST="1.2.18"
VERSION_EGGNOG="5.0.0" # See http://eggnog5.embl.de/#/app/downloads

VERSION_GTDB_COMPLETE= "220.0" # See https://data.gtdb.ecogenomic.org/releases/
VERSION_GTDB_MAIN=VERSION_GTDB_COMPLETE.split('.')[0] # Remove subversion

DB_DIR_PATH = config.get('db_dir')

KEEP_BAM_FILES = is_config_parameter_true(config,'keep_final_coverage_bam_files')

# SAMPLE_NAMES, and POLISH_WITH_SHORT_READS are instantiated in rotary.smk

ANNOTATION_MAP = AnnotationMap(config.get('annotations'))

DFAST_DB_PATH = os.path.join(DB_DIR_PATH,"dfast_" + VERSION_DFAST)

rule download_dfast_db:
    output:
        db_dir=directory(DFAST_DB_PATH),
        install_finished=os.path.join(DB_DIR_PATH,"checkpoints","dfast_" + VERSION_DFAST)
    conda:
        "../envs/annotation_dfast.yaml"
    log:
        "logs/download/dfast_db_download.log"
    benchmark:
        "benchmarks/download/dfast_db_download.txt"
    shell:
        """
        mkdir -p {output.db_dir}
        dfast_file_downloader.py --protein dfast --dbroot {output.db_dir} > {log} 2>&1
        dfast_file_downloader.py --cdd Cog --hmm TIGR --dbroot {output.db_dir} >> {log} 2>&1
        touch {output.install_finished}
        """


rule download_eggnog_db:
    output:
        db=directory(os.path.join(DB_DIR_PATH,"eggnog_" + VERSION_EGGNOG)),
        install_finished=os.path.join(DB_DIR_PATH,"checkpoints","eggnog_" + VERSION_EGGNOG)
    conda:
        "../envs/eggnog.yaml"
    log:
        "logs/download/eggnog_db_download.log"
    benchmark:
        "benchmarks/download/eggnog_db_download.txt"
    params:
        db_dir=os.path.join(DB_DIR_PATH,"eggnog_" + VERSION_EGGNOG)
    shell:
        """
        mkdir -p {params.db_dir}
        download_eggnog_data.py -y -M --data_dir {params.db_dir} > {log} 2>&1
        touch {output.install_finished}
        """


# TODO - if there is an error during download, the initial_download_dir is not deleted during cleanup
rule download_gtdb_db:
    output:
        db=directory(os.path.join(DB_DIR_PATH,"GTDB_" + VERSION_GTDB_COMPLETE)),
        install_finished=os.path.join(DB_DIR_PATH,"checkpoints","GTDB_" + VERSION_GTDB_COMPLETE + "_download")
    log:
        "logs/download/gtdb_db_download.log"
    benchmark:
        "benchmarks/download/gtdb_db_download.txt"
    params:
        db_dir_root=os.path.join(DB_DIR_PATH),
        initial_download_dir=os.path.join(DB_DIR_PATH,"release" + VERSION_GTDB_MAIN),
        db_dir=os.path.join(DB_DIR_PATH,"GTDB_" + VERSION_GTDB_COMPLETE),
        url=f"https://data.ace.uq.edu.au/public/gtdb/data/releases/release{VERSION_GTDB_MAIN}/{VERSION_GTDB_COMPLETE}/"
            f"auxillary_files/gtdbtk_package/full_package/gtdbtk_r{VERSION_GTDB_MAIN}_data.tar.gz"
    shell:
        """
        mkdir -p {params.db_dir_root}

        if [[ -d {params.initial_download_dir} ]]; then
          echo "GTDB cannot be downloaded because of pre-existing dir: {params.initial_download_dir}"
          exit 1
        fi

        wget -O - {params.url} 2> {log} | tar -C {params.db_dir_root} -xzf - >> {log} 2>&1
        mv {params.initial_download_dir} {params.db_dir}

        touch {output.install_finished}
        """


rule validate_gtdb:
    input:
        os.path.join(DB_DIR_PATH,"checkpoints","GTDB_" + VERSION_GTDB_COMPLETE + "_download")
    output:
        os.path.join(DB_DIR_PATH,"checkpoints","GTDB_" + VERSION_GTDB_COMPLETE + "_validate")
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/download/gtdb_db_validate.log"
    benchmark:
        "benchmarks/download/gtdb_db_validate.txt"
    params:
        db_dir=os.path.join(DB_DIR_PATH,"GTDB_" + VERSION_GTDB_COMPLETE)
    shell:
        """
        # Split the "progress bar" style output into multiple lines with sed
        # See: https://stackoverflow.com/a/60786606 (accessed 2022.08.02)
        GTDBTK_DATA_PATH={params.db_dir} gtdbtk check_install | sed 's/\\r/\\n/g' > {log} 2>&1

        touch {output}
        """

rule build_gtdb_mash_ref_database:
    input:
        os.path.join(DB_DIR_PATH,"checkpoints","GTDB_" + VERSION_GTDB_COMPLETE + "_validate")
    output:
        ref_msh_file=os.path.join(DB_DIR_PATH,"GTDB_" + VERSION_GTDB_COMPLETE + '_mash','gtdb_ref_sketch.msh'),
        ref_genome_path_list=temp(os.path.join(DB_DIR_PATH,"GTDB_" + VERSION_GTDB_COMPLETE + '_mash','ref_mash_genomes.txt'))
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/download/gtdb_db_mash.log"
    benchmark:
        "benchmarks/download/gtdb_mash.txt"
    threads:
        config.get("threads",1)
    params:
        fast_ani_genomes_dir=os.path.join(DB_DIR_PATH,"GTDB_" + VERSION_GTDB_COMPLETE,'skani','database')
    shell:
        """
        find {params.fast_ani_genomes_dir} -name *_genomic.fna.gz -type f > {output.ref_genome_path_list}
        mash sketch -l {output.ref_genome_path_list} -p {threads} -o {output.ref_msh_file} -k 16 -s 5000 > {log} 2>&1   
        """

rule download_checkm_db:
    output:
        db=directory(os.path.join(DB_DIR_PATH,"checkm2")),
        install_finished=os.path.join(DB_DIR_PATH,"checkpoints","checkm2")
    conda:
        "../envs/checkm2.yaml"
    log:
        "logs/download/checkm2_db_download.log"
    benchmark:
        "benchmarks/download/checkm2_db_download.txt"
    shell:
        """
        mkdir -p {output.db}
        checkm2 database --download --path {output.db} > {log} 2>&1
        touch {output.install_finished}
        """

rule annotation_download:
    input:
        checkm_install_finished = os.path.join(DB_DIR_PATH,"checkpoints","checkm2"),
        gtdb_ref_msh_file=os.path.join(DB_DIR_PATH,"GTDB_" + VERSION_GTDB_COMPLETE + '_mash','gtdb_ref_sketch.msh'),
        eggnog_install_finished=os.path.join(DB_DIR_PATH,"checkpoints","eggnog_" + VERSION_EGGNOG),
        dfast_install_finished=os.path.join(DB_DIR_PATH,"checkpoints","dfast_" + VERSION_DFAST)
    output:
        annotation_downloaded=touch(os.path.join(DB_DIR_PATH,"checkpoints","annotation_downloaded"))

# TODO - I might need to remove special characters from strain name to use for locus tag prefix
# TODO - can I auto-predict genome completeness, names, types, topologies?
rule run_dfast:
    input:
        contigs="{sample}/circularize/{sample}_circularize.fasta",
        install_finished=os.path.join(DB_DIR_PATH,"checkpoints","dfast_" + VERSION_DFAST) if ANNOTATION_MAP.dfast_func else []
    output:
        dfast_genome="{sample}/annotation/dfast/{sample}_genome.fna",
        dfast_cds="{sample}/annotation/dfast/{sample}_cds.fna",
        dfast_proteins="{sample}/annotation/dfast/{sample}_protein.faa",
        dfast_embl="{sample}/annotation/dfast/{sample}_genome.embl",
        dfast_genbank="{sample}/annotation/dfast/{sample}_genome.gbk",
        dfast_gff="{sample}/annotation/dfast/{sample}_genome.gff",
        dfast_rna="{sample}/annotation/dfast/{sample}_rna.fna",
        dfast_pseudogene="{sample}/annotation/dfast/{sample}_pseudogene_summary.tsv",
        outdir=directory("{sample}/annotation/dfast")
    conda:
        "../envs/annotation_dfast.yaml"
    log:
        "{sample}/logs/annotation/annotation_dfast.log"
    benchmark:
        "{sample}/benchmarks/annotation/annotation_dfast.txt"
    params:
        strain='{sample}',
        no_functional_annotation="--no_func_anno" if not ANNOTATION_MAP.dfast_func else '',
        dfast_db_root = f"--dbroot {DFAST_DB_PATH}" if ANNOTATION_MAP.dfast_func else ''
    threads:
        config.get("threads",1)
    shell:
        """
        dfast --force \
          {params.dfast_db_root} \
          -g {input.contigs} \
          -o {output.outdir} \
          --strain {params.strain} \
          --locus_tag_prefix {params.strain} \
          {params.no_functional_annotation} \
          --cpu {threads} > {log} 2>&1
         # --complete t \
         # --seq_names "Chromosome,unnamed" \
         # --seq_types "chromosome,plasmid" \
         # --seq_topologies "circular,circular" 

         mv {output.outdir}/genome.fna {output.dfast_genome}
         mv {output.outdir}/cds.fna {output.dfast_cds}
         mv {output.outdir}/protein.faa {output.dfast_proteins}
         mv {output.outdir}/genome.embl {output.dfast_embl}
         mv {output.outdir}/genome.gbk {output.dfast_genbank}
         mv {output.outdir}/genome.gff {output.dfast_gff}
         mv {output.outdir}/rna.fna {output.dfast_rna}
         mv {output.outdir}/pseudogene_summary.tsv {output.dfast_pseudogene}
         """


# TODO - add option to control whether --dbmem flag is set (uses more RAM but does faster analysis)
rule run_eggnog:
    input:
        protein="{sample}/annotation/dfast/{sample}_protein.faa",
        install_finished=os.path.join(DB_DIR_PATH,"checkpoints","eggnog_" + VERSION_EGGNOG)
    output:
        eggnog_annotations="{sample}/annotation/eggnog/{sample}.emapper.annotations",
        outdir=directory("{sample}/annotation/eggnog")
    conda:
        "../envs/eggnog.yaml"
    log:
        "{sample}/logs/annotation/eggnog.log"
    benchmark:
        "{sample}/benchmarks/annotation/eggnog.txt"
    params:
        prefix='{sample}',
        db=directory(os.path.join(DB_DIR_PATH,"eggnog_" + VERSION_EGGNOG)),
        sensmode=config.get("eggnog_sensmode"),
        search_tool=config.get('eggnog_search_tool')
    threads:
        config.get("threads",1)
    shell:
        """
        mkdir -p {output.outdir}/tmp
        emapper.py --cpu {threads} -i {input.protein} --itype proteins -m {params.search_tool} \
          --sensmode {params.sensmode} --dbmem --output {params.prefix} --output_dir {output.outdir} \
          --temp_dir {output.outdir}/tmp \
          --data_dir {params.db} --override > {log} 2>&1
        rm -r {output.outdir}/tmp
        """


rule run_gtdbtk:
    input:
        genome="{sample}/annotation/dfast/{sample}_genome.fna",
        setup_finished=os.path.join(DB_DIR_PATH,"checkpoints","GTDB_" + VERSION_GTDB_COMPLETE + "_validate"),
        ref_msh_file=os.path.join(DB_DIR_PATH,"GTDB_" + VERSION_GTDB_COMPLETE + '_mash','gtdb_ref_sketch.msh')
    output:
        batchfile=temp("{sample}/annotation/gtdbtk/batchfile.tsv"),
        annotation="{sample}/annotation/gtdbtk/{sample}_gtdbtk.summary.tsv",
        outdir=directory("{sample}/annotation/gtdbtk/run_files")
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "{sample}/logs/annotation/gtdbtk.log"
    benchmark:
        "{sample}/benchmarks/annotation/gtdbtk.txt"
    params:
        db_dir=directory(os.path.join(DB_DIR_PATH,"GTDB_" + VERSION_GTDB_COMPLETE)),
        genome_id='{sample}',# Get sample name from wildcards.
        gtdbtk_mode="--full_tree" if config.get("gtdbtk_mode") == "full_tree" else ""
    threads:
        config.get("threads",1)
    shell:
        """
        printf "{input.genome}\t{params.genome_id}\n" > {output.batchfile}
        GTDBTK_DATA_PATH={params.db_dir} \
        gtdbtk classify_wf --batchfile {output.batchfile} --out_dir {output.outdir} {params.gtdbtk_mode} \
          --mash_db {input.ref_msh_file} --cpus {threads} --pplacer_cpus {threads} > {log} 2>&1
        head -n 1 {output.outdir}/gtdbtk.*.summary.tsv | sort -u > {output.annotation}
        tail -n +2 {output.outdir}/gtdbtk.*.summary.tsv >> {output.annotation}
        """

rule aggregate_gtdbtk_reports:
    input:
        expand("{sample}/annotation/gtdbtk/{sample}_gtdbtk.summary.tsv", sample=SAMPLE_NAMES)
    output:
        'aggregate_stats/annotation/aggregate_gtdbtk_summary_report.tsv'
    benchmark:
        "benchmarks/annotation/aggregate_gtdbtk.benchmark.txt"
    run:
        combine_tabular_reports(input,output[0])


rule run_checkm2:
    input:
        genome="{sample}/annotation/dfast/{sample}_genome.fna",
        install_finished=os.path.join(DB_DIR_PATH,"checkpoints","checkm2")
    output:
        quality_report="{sample}/annotation/checkm/{sample}_checkm_quality_report.tsv",
        outdir=directory("{sample}/annotation/checkm/")
    log:
        "{sample}/logs/annotation/checkm.log"
    benchmark:
        "{sample}/benchmarks/annotation/checkm.txt"
    conda:
        "../envs/checkm2.yaml"
    params:
        db=directory(os.path.join(DB_DIR_PATH,"checkm2","CheckM2_database","uniref100.KO.1.dmnd"))
    threads:
        config.get("threads",1)
    shell:
        """
        mkdir -p {output.outdir}
        checkm2 predict --threads {threads} --input {input.genome} --output-directory {output.outdir} \
            --database_path {params.db} > {log} 2>&1
        mv {output.outdir}/quality_report.tsv {output.quality_report}
        """

rule aggregate_checkm2_reports:
    input:
        expand("{sample}/annotation/checkm/{sample}_checkm_quality_report.tsv", sample=SAMPLE_NAMES)
    output:
        'aggregate_stats/annotation/aggregate_checkm_quality_report.tsv'
    benchmark:
        "benchmarks/annotation/aggregate_checkm.benchmark.txt"
    run:
        combine_tabular_reports(input,output[0])


if POLISH_WITH_SHORT_READS == True:
    # TODO - clarify name compared to previous mapping step
    rule calculate_final_short_read_coverage:
        input:
            qc_short_r1="{sample}/qc/{sample}_qc_R1.fastq.gz",
            qc_short_r2="{sample}/qc/{sample}_qc_R2.fastq.gz",
            dfast_genome="{sample}/annotation/dfast/{sample}_genome.fna"
        output:
            mapping="{sample}/annotation/coverage/{sample}_short_read.bam" if KEEP_BAM_FILES else temp("{sample}/annotation/coverage/{sample}_short_read.bam"),
            index=temp("{sample}/annotation/coverage/{sample}_short_read.bam.bai"),
            coverage=temp("{sample}/annotation/coverage/{sample}_short_read_coverage_raw.tsv"),
            read_mapping_files= temp(multiext("{sample}/annotation/dfast/{sample}_genome.fna",
                *READ_MAPPING_FILE_EXTENSIONS)) # Variable declared in polish.smk
        conda:
            "../envs/mapping.yaml"
        log:
            "{sample}/logs/annotation/calculate_final_short_read_coverage.log"
        benchmark:
            "{sample}/benchmarks/annotation/calculate_final_short_read_coverage.txt"
        params:
            mem_per_thread = int(config.get("memory") / config.get("threads",1))
        threads:
            config.get("threads",1)
        resources:
            mem=config.get("memory")
        shell:
            """
            bwa-mem2 index {input.dfast_genome} 2> {log}
            bwa-mem2 mem -t {threads} {input.dfast_genome} {input.qc_short_r1} {input.qc_short_r2} 2>> {log} | \
              samtools view -b -@ {threads} 2>> {log} | \
              samtools sort -@ {threads} -m {params.mem_per_thread}G 2>> {log} \
              > {output.mapping}
            samtools index -@ {threads} {output.mapping}
            samtools coverage {output.mapping} > {output.coverage}
            """


rule calculate_final_long_read_coverage:
    input:
        contigs="{sample}/annotation/dfast/{sample}_genome.fna",
        qc_long_reads="{sample}/qc/{sample}_qc_long.fastq.gz"
    output:
        mapping="{sample}/annotation/coverage/{sample}_long_read.bam" if KEEP_BAM_FILES else temp("{sample}/annotation/coverage/{sample}_long_read.bam"),
        index=temp("{sample}/annotation/coverage/{sample}_long_read.bam.bai"),
        coverage=temp("{sample}/annotation/coverage/{sample}_long_read_coverage_raw.tsv")
    conda:
        "../envs/mapping.yaml"
    log:
        "{sample}/logs/annotation/calculate_final_long_read_coverage.log"
    benchmark:
        "{sample}/benchmarks/annotation/calculate_final_long_read_coverage.txt"
    params:
        mem_per_thread = int(config.get("memory") / config.get("threads",1))
    threads:
        config.get("threads",1)
    resources:
        mem=config.get("memory")
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.contigs} {input.qc_long_reads} 2> {log} | \
          samtools view -b -@ {threads} 2>> {log} | \
          samtools sort -@ {threads} -m {params.mem_per_thread}G 2>> {log} \
          > {output.mapping}
        samtools index -@ {threads} {output.mapping}
        samtools coverage {output.mapping} > {output.coverage}
        """


rule reformat_final_read_coverage:
    input:
        "{sample}/annotation/coverage/{sample}_{long_or_short}_read_coverage_raw.tsv"
    output:
        "{sample}/annotation/coverage/{sample}_{long_or_short}_read_coverage.tsv"
    params:
        sample_id="{sample}",
        read_type="{long_or_short}"
    run:
        read_coverage = pd.read_csv(input[0], sep='\t')\
            .rename(columns={'#rname':'contig'})

        read_coverage.insert(0, column='sample', value=params.sample_id)
        read_coverage.insert(1, column='read_type', value=params.read_type)

        read_coverage.to_csv(output[0], sep='\t', index=False)


rule aggregate_final_read_coverage:
    input:
        expand("{sample}/annotation/coverage/{sample}_{long_or_short}_read_coverage.tsv", sample=SAMPLE_NAMES, long_or_short=['long','short']) \
            if POLISH_WITH_SHORT_READS == True else expand("{sample}/annotation/coverage/{sample}_long_read_coverage.tsv", sample=SAMPLE_NAMES)
    output:
        'aggregate_stats/annotation/aggregate_read_coverage.tsv'
    benchmark:
        "benchmarks/annotation/aggregate_read_coverage.txt"
    run:
        combine_tabular_reports(input,output[0])


rule symlink_logs:
    input:
        dfast_genome="{sample}/annotation/dfast/{sample}_genome.fna"
    output:
        logs=temp(directory("{sample}/annotation/logs")),
        stats=temp(directory("{sample}/annotation/stats"))
    params:
        logs="{sample}/logs",
        stats="{sample}/stats"
    run:
        source_relpath = os.path.relpath(str(params.logs),os.path.dirname(str(output.logs)))
        os.symlink(source_relpath,str(output.logs))

        source_relpath = os.path.relpath(str(params.stats),os.path.dirname(str(output.stats)))
        os.symlink(source_relpath,str(output.stats))


rule summarize_annotation:
    input:
        "{sample}/annotation/dfast/{sample}_genome.fna",
        "{sample}/annotation/eggnog/{sample}.emapper.annotations" if ANNOTATION_MAP.eggnog else [],
        "{sample}/annotation/gtdbtk/{sample}_gtdbtk.summary.tsv" if ANNOTATION_MAP.gtdbtk else [],
        "{sample}/annotation/checkm/" if ANNOTATION_MAP.checkm2 else [],
        "{sample}/annotation/coverage/{sample}_short_read_coverage.tsv" if POLISH_WITH_SHORT_READS and ANNOTATION_MAP.coverage else [],
        "{sample}/annotation/coverage/{sample}_long_read_coverage.tsv" if ANNOTATION_MAP.coverage else [],
        "{sample}/annotation/logs",
        "{sample}/annotation/stats"
    output:
        "{sample}/{sample}_annotation_summary.zip"
    log:
        "{sample}/logs/annotation/summarize_annotation.log"
    params:
        zipdir="{sample}/annotation"
    shell:
        """
        zip -r {output} {params.zipdir}/* -x \*.bam\* gtdbtk/run_files/\* > "{log}" 2>&1
        """


rule annotation:
    input:
        summeries=expand("{sample}/{sample}_annotation_summary.zip",sample=SAMPLE_NAMES),
        aggregate_checkm_report='aggregate_stats/annotation/aggregate_checkm_quality_report.tsv',
        aggregate_gtdbtk_report='aggregate_stats/annotation/aggregate_gtdbtk_summary_report.tsv',
        aggregate_final_read_coverage='aggregate_stats/annotation/aggregate_read_coverage.tsv'
    output:
        temp(touch("checkpoints/annotation"))
