import pandas as pd
import seaborn as sns
from scipy.stats import kruskal
import matplotlib.pyplot as plt
import numpy as np
from statsmodels.stats.multitest import multipletests
import streamlit as st

def read_input_files(gene_matrix_file, metadata_file, annotation_file):
    gene_matrix = pd.read_csv(gene_matrix_file, index_col=0, sep="\t")
    metadata = pd.read_csv(metadata_file, index_col=0, sep=None, engine='python')
    annotation = pd.read_csv(annotation_file, index_col=0, sep="\t")
    return gene_matrix, metadata, annotation

def differential_abundance_test(gene_matrix, metadata, annotation, metadata_column, annotation_column, threshold, lfc_threshold):
    # Merge metadata and gene_matrix on the specified metadata column
    merged_data = pd.merge(annotation[[annotation_column]], gene_matrix, left_index=True, right_index=True)

    conditions = (merged_data[annotation_column].notna()) & (merged_data[annotation_column] != '') & (merged_data[annotation_column] != ' ') & (merged_data[annotation_column] != '-') & (merged_data[annotation_column] != '_')

    # Apply the filter to keep rows that meet the conditions
    merged_data = merged_data[conditions]

    merged_data[annotation_column] = merged_data[annotation_column].str.split(',')

    # Explode the rows to duplicate entries when there are multiple values
    merged_data = merged_data.explode(annotation_column)

    # Collapse gene matrix by annotation column
    collapsed_data = merged_data.groupby(annotation_column).sum()

    # Transpose the collapsed data to make samples as rows
    transposed_data = collapsed_data.T

    # Merge metadata with transposed_data based on sample names
    merged_metadata = pd.merge(metadata[[metadata_column]], transposed_data, left_index=True, right_index=True)

    conditions = (merged_metadata[metadata_column].notna()) & (merged_metadata[metadata_column] != '') & (merged_metadata[metadata_column] != ' ') & (merged_metadata[metadata_column] != '-') & (merged_metadata[metadata_column] != '_')

    # Apply the filter to keep rows that meet the conditions
    merged_metadata = merged_metadata[conditions]

    # Group by the specified metadata column
    grouped_data_mean = merged_metadata.groupby(metadata_column).mean().T

    # Add a line to calculate the log-transformed mean values
    grouped_data_log_mean = np.log1p(grouped_data_mean)

    # Concatenate both dataframes
    grouped_data = pd.concat([grouped_data_mean, grouped_data_log_mean.add_suffix('_log')], axis=1)

    # Calculate max log fold change 
    grouped_data['MaxLog2FC'] = grouped_data_log_mean.apply(lambda row: max(row) - min(row), axis=1)

    abundance_columns = grouped_data.index

    groups = merged_metadata[metadata_column].unique()
    
    for abundance_column in abundance_columns:
        abundance_groups = [merged_metadata[abundance_column][merged_metadata[metadata_column] == group] for group in groups]

        try:
            statistic, p_value = kruskal(*abundance_groups)
        except ValueError as e:
            print(f"Warning: {e}. Setting both statistic and p_value to 1.")
            statistic = 1
            p_value = 1

        # Access grouped_data as a DataFrame to assign the p-value
        grouped_data.at[abundance_column, 'P-value'] = p_value

    _, corrected_p_values, _, _ = multipletests(grouped_data['P-value'], method='fdr_bh')

    # Add corrected p-values to the results dictionary
    grouped_data['Corrected P-value'] = corrected_p_values

    # Save significant gene information to a file
    output_file = f'abundance_{metadata_column}_{annotation_column}.csv'
    grouped_data.to_csv(output_file, index=True)

    # Filter by LFC and p-value
    filtered_data = grouped_data[(grouped_data['MaxLog2FC'] >= lfc_threshold) & (grouped_data['Corrected P-value'] <= threshold)]

    # Return only the log abundance for each group column
    result_matrix = filtered_data[groups]

    return result_matrix

def plot_heatmap(result_matrix, metadata_column, annotation_column):
    # Calculate the size based on the number of rows and columns
    height, width = result_matrix.shape[0] * 0.2 + 4, result_matrix.shape[1] * 0.2 + 1.5

    # Define the color palette with continuous blending
    cmap = sns.color_palette("crest", as_cmap=True)

    # Create the clustermap with xticklabels rotated by 90 degrees
    g = sns.clustermap(np.log2(result_matrix+0.01), col_cluster=True, row_cluster=True, xticklabels=True, figsize=(width, height), cmap=cmap, fmt=".2f", linewidths=.5, )

    # Rotate y-axis labels to ensure they are all shown
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)

    # Increase the resolution for high-quality output
    g.savefig(f'heatmap_{metadata_column}_{annotation_column}.png', bbox_inches='tight', dpi=300)

    # Display heatmap in Streamlit app
    st.pyplot()

def abundance():
    st.title("Differential Abundance Analysis")

    # File Uploads
    st.sidebar.header("Upload Files")
    metadata_file = st.sidebar.file_uploader("Select Metadata File (TSV)", type=["tsv", "txt"])
    gene_matrix_file = st.sidebar.file_uploader("Select Gene Matrix File (TSV)", type=["tsv", "txt"])
    annotation_file = st.sidebar.file_uploader("Select Annotation File (TSV)", type=["tsv", "txt"])

    # Metadata and Annotation Columns
    st.sidebar.header("Select Metadata and Annotation Columns")
    metadata_columns = []
    annotation_columns = []

    if gene_matrix_file and metadata_file and annotation_file:
        gene_matrix, metadata, annotation = read_input_files(gene_matrix_file, metadata_file, annotation_file)
        metadata_columns = list(metadata.columns)
        annotation_columns = list(annotation.columns)

    metadata_column = st.sidebar.selectbox("Select Metadata Column", metadata_columns)
    annotation_column = st.sidebar.selectbox("Select Annotation Column", annotation_columns)

    # Thresholds
    st.sidebar.header("Set Thresholds")
    pval = st.sidebar.slider("Select Significance Threshold (P-value)", min_value=0.0001, max_value=0.1, value=0.01, step=0.01)
    lfc = st.sidebar.slider("Select Log Fold Change Threshold", min_value=0.1, max_value=2.0, value=1.0, step=0.1)

    # Run Analysis
    if st.button("Run Analysis"):
        if gene_matrix_file and metadata_file and annotation_file:
            result_matrix = differential_abundance_test(gene_matrix, metadata, annotation, metadata_column, annotation_column, pval, lfc)
            plot_heatmap(result_matrix, metadata_column, annotation_column)
            st.success("Analysis completed successfully! Check the output files.")

if __name__ == "__main__":
    abundance()