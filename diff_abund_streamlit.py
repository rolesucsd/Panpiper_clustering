import pandas as pd
import seaborn as sns
from scipy.stats import kruskal
import matplotlib.pyplot as plt
import numpy as np
from statsmodels.stats.multitest import multipletests
import streamlit as st

def read_input_files(gene_matrix_file, metadata_file, annotation_file):
    try:
        gene_matrix = pd.read_csv(gene_matrix_file, index_col=0, sep="\t")
        metadata = pd.read_csv(metadata_file, index_col=0, sep=None, engine='python')
        annotation = pd.read_csv(annotation_file, index_col=0, sep="\t")
        return gene_matrix, metadata, annotation
    except Exception as e:
        st.error(f"Error reading input files: {str(e)}")
        return None, None, None

def differential_abundance_test(gene_matrix, metadata, annotation, metadata_column, annotation_column, threshold, lfc_threshold):
    try:
        # ... (your existing code)

        return result_matrix
    except Exception as e:
        st.error(f"Error in differential abundance test: {str(e)}")
        return None

def plot_heatmap(result_matrix, metadata_column, annotation_column):
    try:
        # ... (your existing code)

        st.pyplot()
    except Exception as e:
        st.error(f"Error plotting heatmap: {str(e)}")

def abundance():
    try:
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
            if gene_matrix is not None and metadata is not None and annotation is not None:
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
                if result_matrix is not None:
                    plot_heatmap(result_matrix, metadata_column, annotation_column)
                    st.success("Analysis completed successfully! Check the output files.")
    except Exception as e:
        st.error(f"An unexpected error occurred: {str(e)}")

if __name__ == "__main__":
    abundance()
