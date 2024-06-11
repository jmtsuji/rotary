#!/usr/bin/env python

"""
Created by: Lee Bergstrand (2024)

Description: Code for generating QC reports.
"""

from zipfile import ZipFile

import numpy as np
import pandas as pd


def extract_table_data_from_multiqc_zip(multiqc_zip_path, table_filename='multiqc_fastqc.txt', delimiter='\t'):
    """
    Extract table data from a multiqc zip file.

    :param multiqc_zip_path: Path to the multiqc zip file.
    :param table_filename: Name of the table file inside the multiqc zip file. Default is 'multiqc_fastqc.txt'.
    :param delimiter: Delimiter used in the table file. Default is '\t'.

    :return: A pandas DataFrame containing the table file data from the multiqc zip file.
    """
    with ZipFile(multiqc_zip_path) as zip_ref:
        for filename in zip_ref.namelist():
            if filename == table_filename:
                with zip_ref.open(filename) as fastqc_data:
                    dataframe = pd.read_csv(fastqc_data, delimiter=delimiter)
                return dataframe


def sanatize_fastqc_dataframe(raw_fastqc_dataframe):
    raw_sample_data_column_name = 'Raw_Sample_Data'

    # Rename the sample column to a more general name because it will be split into multiple columns.
    raw_fastqc_dataframe.rename(columns={'Sample': raw_sample_data_column_name}, inplace=True)

    # The QC file type designators that are included in the QC file names.
    file_types = ['_adapter_trim', '_reformat', '_quality_trim', '_reformat', '_filter']

    generate_fastqc_sample_name_stage_and_direction_columns(raw_fastqc_dataframe, raw_sample_data_column_name,
                                                            file_types)

    convert_string_fastqc_data_to_numeric(raw_fastqc_dataframe)

    cleaned_fastqc_dataframe = raw_fastqc_dataframe[['Sample', 'Stage', 'Direction', 'Total Sequences',
                                                     'Total Bases (Mbp)', 'Sequence Length Min', 'Sequence Length Max',
                                                     'avg_sequence_length', 'median_sequence_length', '%GC']]

    return cleaned_fastqc_dataframe


def generate_fastqc_sample_name_stage_and_direction_columns(fastqc_dataframe, sample_column_name, file_types):
    sample_name_no_direction_column_name = 'Sample_Name_No_Direction'

    # Split the sample_name column so that the direction of the sample is in its own column.
    processed_sample_columns = fastqc_dataframe[sample_column_name].str.rsplit('_', n=1, expand=True)
    processed_sample_columns.columns = [sample_name_no_direction_column_name, 'Direction']
    sample_column_no_direction = processed_sample_columns[sample_name_no_direction_column_name]

    # Remove QC file type information from the sample names.
    for qc_file_type_designators in file_types:
        sample_column_no_direction = sample_column_no_direction.str.replace(qc_file_type_designators, '', regex=False)
    processed_sample_columns['Sample'] = sample_column_no_direction

    # Generate the QC stage of each sample by removing the sample name
    # from the recently generated sample_column_no_direction column.
    sample_types = processed_sample_columns.apply(
        lambda row: row[sample_name_no_direction_column_name].replace(row['Sample'], ''), axis=1)
    processed_sample_columns['Stage'] = sample_types.replace('', 'raw')  # Samples with no sample type should be raw.
    processed_sample_columns.drop(columns=sample_name_no_direction_column_name, inplace=True)
    fastqc_dataframe[['Direction', 'Sample', 'Stage']] = processed_sample_columns


def convert_string_fastqc_data_to_numeric(fastqc_dataframe):
    # Create a new total bases column that uses float values instead of strings.
    fastqc_dataframe['Total Bases (Mbp)'] = fastqc_dataframe['Total Bases'].str.replace(' Mbp', ' ').astype(
        np.float64)

    # Split the sequence length column into two columns that two float columns.
    fastqc_dataframe[['Sequence Length Min', 'Sequence Length Max']] = fastqc_dataframe['Sequence length'].str.split(
        '-', expand=True).astype(np.float64)


def generate_fastqc_comparison_table(raw_multiqc_dataframe, direction, start_qc_file_type, end_qc_file_type):
    cleaned_fastqc_dataframe = sanatize_fastqc_dataframe(raw_multiqc_dataframe)

    direction_column_name = 'Direction'
    qc_stage_column_name = 'Stage'
    name_column_name = 'Sample'
    before_qc_name = 'Before QC'
    after_qc_name = 'After QC'

    # Filter down the comparison table to a certain FASTQ file direction.
    single_direction_table = cleaned_fastqc_dataframe[cleaned_fastqc_dataframe[direction_column_name] == direction]
    single_direction_table = single_direction_table.drop(columns=direction_column_name)

    # Filter down the comparison table to just the before and after qc file types.
    # Replace these file type names with before and after QC.
    comparison_table = single_direction_table[single_direction_table[qc_stage_column_name].isin(
        [start_qc_file_type, end_qc_file_type])].replace(
        {start_qc_file_type: before_qc_name, end_qc_file_type: after_qc_name})

    before_qc_samples = comparison_table[comparison_table[qc_stage_column_name] == before_qc_name].set_index(
        name_column_name)
    after_qc_samples = comparison_table[comparison_table[qc_stage_column_name] == after_qc_name].set_index(
        name_column_name)

    before_qc_samples.drop(columns=qc_stage_column_name, inplace=True)
    after_qc_samples.drop(columns=qc_stage_column_name, inplace=True)

    percentage_change = (after_qc_samples - before_qc_samples) / before_qc_samples * 100
    percentage_change[qc_stage_column_name] = 'Change (%)'

    comparison_table.set_index(name_column_name, inplace=True)
    final_data = pd.concat([comparison_table, percentage_change]).sort_index().reset_index().set_index(
        [name_column_name, qc_stage_column_name]).transpose()

    return final_data
