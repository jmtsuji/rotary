# rotary: utilities and workflow for long-read DNA assemblies including circular elements
# Copyright Jackson M. Tsuji and Lee H. Bergstrand, 2023

import os
from snakemake.utils import min_version

from pungi.dataset import generate_dataset_from_sample_tsv
from pungi.utils import symlink_or_compress, is_config_parameter_true

SAMPLE_TSV_PATH = 'samples.tsv'

# When running rotary in download mode, there can be situations where samples.tsv is not present in the CWD.
# We bypass processing samples if no samples.tsv file is found in the CWD.
if os.path.isfile(os.path.join(os.getcwd(), SAMPLE_TSV_PATH)):
    SAMPLES = generate_dataset_from_sample_tsv(SAMPLE_TSV_PATH)

    if SAMPLES.files_per_sample==3:
        DATASET_HAS_SHORT_READS=True
    else:
        DATASET_HAS_SHORT_READS=False

    SAMPLE_NAMES = list(SAMPLES.identifiers)
else:
    # Some variables are required for evaluating the DAG because they are used in rules.
    # We mock them if no samples.tsv file is found.
    SAMPLE_NAMES=['mock_sample_name']
    DATASET_HAS_SHORT_READS=False

# Specify the minimum snakemake version allowable
min_version("7.0")
# Specify shell parameters
shell.executable("/bin/bash")
shell.prefix("set -o pipefail; ")

POLISH_WITH_SHORT_READS = is_config_parameter_true(config,'polish_with_short_reads')

if DATASET_HAS_SHORT_READS is False:
    POLISH_WITH_SHORT_READS=False

rule all:
    input:
        "checkpoints/qc",
        "checkpoints/assembly",
        "checkpoints/polish",
        "checkpoints/circularize",
        "checkpoints/annotation"

rule set_up_sample_directories:
    input:
        SAMPLE_TSV_PATH
    output:
        long_reads = expand("{sample}/raw/{sample}_long.fastq.gz", sample=SAMPLE_NAMES),
        short_R1_reads = expand("{sample}/raw/{sample}_R1.fastq.gz", sample=SAMPLE_NAMES) if DATASET_HAS_SHORT_READS else [],
        short_R2_reads = expand("{sample}/raw/{sample}_R2.fastq.gz", sample=SAMPLE_NAMES) if DATASET_HAS_SHORT_READS else []
    run:
        for sample in SAMPLES:
            identifier = sample.identifier
            symlink_or_compress(sample.long_read_path,f'{identifier}/raw/{identifier}_long.fastq.gz')

            if DATASET_HAS_SHORT_READS:
                symlink_or_compress(sample.short_read_left_path,f'{identifier}/raw/{identifier}_R1.fastq.gz')
                symlink_or_compress(sample.short_read_right_path,f'{identifier}/raw/{identifier}_R2.fastq.gz')

# Include various modules.
include: './qc.smk'
include: './assembly.smk'
include: './polish.smk'
include: './circularize.smk'
include: './annotation.smk'

DB_DIR_PATH = config.get('db_dir')

rule download:
    input:
        annotation_downloaded=os.path.join(DB_DIR_PATH,"checkpoints","annotation_downloaded"),
        circularize_downloaded=os.path.join(DB_DIR_PATH,"checkpoints","circularize_downloaded"),
        qc_downloaded=os.path.join(DB_DIR_PATH,"checkpoints","qc_downloaded")

# TODO - add nice summaries of reads removed during QC, polishing stats, assembly stats, circular vs. non.
#        These aren't in the final ZIP summary at the moment.
