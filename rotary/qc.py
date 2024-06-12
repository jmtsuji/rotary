#!/usr/bin/env python

"""
Created by: Lee Bergstrand (2024)

Description: Code for generating QC reports.
"""
from zipfile import ZipFile

import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype

PERCENTAGE_CHANGE_NAME = 'Change (%)'
BEFORE_QC_NAME = 'Before QC'
AFTER_QC_NAME = 'After QC'

before_after_categorical_type = CategoricalDtype(categories=[BEFORE_QC_NAME, AFTER_QC_NAME, PERCENTAGE_CHANGE_NAME],
                                                 ordered=True)


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


def sanitize_fastqc_data(raw_fastqc_data: pd.DataFrame, file_types=None, data_columns_to_include=None):
    """
    This method performs various cleaning operations on a raw FastQC DataFrame.
    Various columns are converted from strings to numeric.
    The sample column is split into three columns for sample name, stage, and direction.

    :param raw_fastqc_data: The raw FastQC DataFrame that needs to be sanitized.
    :param file_types: The strings used in the file names to designate different QC file types.
    :param data_columns_to_include: The data columns to include in sanitized DataFrame.
    :return: The cleaned FastQC DataFram that contains only the desired columns.
    """
    if data_columns_to_include is None:
        data_columns_to_include = ['Sample', 'Stage', 'Direction', 'Total Sequences',
                                   'Total Bases (Mbp)', 'Sequence Length Min', 'Sequence Length Max',
                                   'avg_sequence_length', 'median_sequence_length', '%GC']
    if file_types is None:
        file_types = ['_adapter_trim', '_reformat', '_quality_trim', '_reformat', '_filter']

    # Rename the sample column to a more general name because it will be split into multiple columns.
    sample_data_column_name = 'Raw_Sample_Data'
    raw_fastqc_data.rename(columns={'Sample': sample_data_column_name}, inplace=True)

    generate_fastqc_sample_name_stage_and_direction_columns(raw_fastqc_data, sample_data_column_name,
                                                            file_types)
    convert_string_fastqc_data_to_numeric(raw_fastqc_data)

    return raw_fastqc_data[data_columns_to_include]


def generate_fastqc_sample_name_stage_and_direction_columns(fastqc_data: pd.DataFrame, sample_column, file_types):
    """
    Splits the input sample column into sample name, type, and direction columns.

    :param fastqc_data: A pandas DataFrame containing the FastQC data.
    :param sample_column: The name of the column in the DataFrame that contains the raw sample names.
    :param file_types: A list of file type strings to remove from the sample names.
    """
    temp_sample_no_direction_col_name = 'Sample_Name_No_Direction'

    # Split the sample_name column so that the direction of the sample is in its own column.
    processed_sample_columns = fastqc_data[sample_column].str.rsplit('_', n=1, expand=True)
    processed_sample_columns.columns = [temp_sample_no_direction_col_name, 'Direction']
    sample_column_no_direction = processed_sample_columns[temp_sample_no_direction_col_name]

    # Remove QC file type information from the sample names.
    for qc_file_type_designators in file_types:
        sample_column_no_direction = sample_column_no_direction.str.replace(qc_file_type_designators, '', regex=False)
    processed_sample_columns['Sample'] = sample_column_no_direction

    # Generate the QC stage of each sample by removing the sample name
    # from the recently generated sample_column_no_direction column.
    sample_types = processed_sample_columns.apply(
        lambda row: row[temp_sample_no_direction_col_name].replace(row['Sample'], ''), axis=1)
    processed_sample_columns['Stage'] = sample_types.replace('', 'raw')  # Samples with no sample type should be raw.
    processed_sample_columns.drop(columns=temp_sample_no_direction_col_name, inplace=True)
    fastqc_data[['Direction', 'Sample', 'Stage']] = processed_sample_columns


def convert_string_fastqc_data_to_numeric(fastqc_data: pd.DataFrame, total_bases_column='Total Bases',
                                          sequence_length_column='Sequence length'):
    """
    Converts total bases and sequence length columns of FastQC data from string format to numeric format.
    The total bases column is directly converted to numeric,
    and the sequence length column is split into min and max sequence length columns.

    :param fastqc_data: A pandas DataFrame containing the FastQC data.
    :param total_bases_column: The name of column in the DataFrame containing total bases data.
    :param sequence_length_column: The name of column in the DataFrame containing sequence length data.
    """
    # Create a new total bases column that uses float values instead of strings.
    fastqc_data['Total Bases (Mbp)'] = fastqc_data[total_bases_column].str.replace(' Mbp', ' ').astype(
        np.float64)

    # Split the sequence length column into two columns that two float columns.
    fastqc_data[['Sequence Length Min', 'Sequence Length Max']] = fastqc_data[
        sequence_length_column].str.split('-', expand=True).astype(np.float64)


def generate_fastqc_before_and_after_comparison_data(single_direction_fastqc_data: pd.DataFrame,
                                                     start_file_type, end_file_type):
    """
    Generates a comparison from the FastQC data of sequencing data stats before and after QC was performed.
    Calculates the percentage change of various sequencing data stats before and after QC.

    :param single_direction_fastqc_data: A DataFrame containing FastQC data for a single direction.
    :param start_file_type: A string representing the start QC file type.
    :param end_file_type: A string representing the end QC file type.
    :return: DataFrame containing sequencing data stats before and after QC was performed.
    """
    qc_stage_column_name = 'Stage'
    name_column_name = 'Sample'

    # Filter down the comparison table to just the before and after qc file types.
    # Replace these file type names with before and after QC.
    before_and_after_qc_table = single_direction_fastqc_data[single_direction_fastqc_data[qc_stage_column_name].isin(
        [start_file_type, end_file_type])].replace(
        {start_file_type: BEFORE_QC_NAME, end_file_type: AFTER_QC_NAME}).set_index(name_column_name)

    # Split the single direction DataFrame into before and after QC DataFrames.
    before_qc_samples = before_and_after_qc_table[before_and_after_qc_table[qc_stage_column_name] == BEFORE_QC_NAME]
    after_qc_samples = before_and_after_qc_table[before_and_after_qc_table[qc_stage_column_name] == AFTER_QC_NAME]

    after_qc_numeric = after_qc_samples.drop(columns=qc_stage_column_name)
    before_qc_numeric = before_qc_samples.drop(columns=qc_stage_column_name)

    # Calculate percentage change as a new dataframe.
    percentage_change = (after_qc_numeric - before_qc_numeric) / before_qc_samples * 100
    percentage_change[qc_stage_column_name] = PERCENTAGE_CHANGE_NAME

    # Join before, after, and percentage change data into a single dataframe.
    final_data = pd.concat([before_and_after_qc_table, percentage_change])
    # Apply a new ordered categorical type that goes before, after, then percentage change.
    final_data[qc_stage_column_name] = final_data[qc_stage_column_name].astype(before_after_categorical_type)

    # Add the Stage column to the index, sort the DataFrame and then return it.
    return final_data.set_index(qc_stage_column_name, append=True).sort_index()


def generate_single_direction_fastqc_data(sanitized_fastqc_data: pd.DataFrame, direction):
    """
    Generate a DataFrame containing only the data for a single FASTQ direction from the cleaned FASTQC data.

    :param sanitized_fastqc_data: The cleaned FASTQC DataFrame, containing information for one or more directions.
    :param direction: The direction to filter the data to.
    :return: A DataFrame containing only the data for the specified direction.
    """
    direction_column_name = 'Direction'
    # Filter down the comparison table to a certain FASTQ file direction.
    single_direction_table = sanitized_fastqc_data[sanitized_fastqc_data[direction_column_name] == direction]
    single_direction_table = single_direction_table.drop(columns=direction_column_name)
    return single_direction_table


def write_fastqc_summary_tsv(sanitized_fastqc_data, output_tsv_path, read_direction):
    """
    Writes summary TSV files for FastQC data.

    :param sanitized_fastqc_data: The cleaned FASTQC DataFrame, containing information for one or more directions.
    :param output_tsv_path: The output path to write the TSV to.
    :param read_direction: Type of read ('R1', 'R2', or 'long').
    """

    single_direction_fastqc_data = generate_single_direction_fastqc_data(sanitized_fastqc_data, read_direction)
    comparison_data = generate_fastqc_before_and_after_comparison_data(single_direction_fastqc_data, 'raw', '_filter')
    comparison_data.to_csv(output_tsv_path, sep='\t', index=True)


def extract_and_combine_sample_fastqc_multiqc_data(sample_fastqc_multiqc_paths):
    """
    Extract the FastQC data from multiple MultiQC zip files and combine it into a single DataFrame.

    :param sample_fastqc_multiqc_paths: List of file paths to the MultiQC reports for each sample.
    :return: Combined FastQC DataFrame for multiple samples.
    """
    return pd.concat([extract_table_data_from_multiqc_zip(path) for path in sample_fastqc_multiqc_paths])



