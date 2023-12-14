import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def set_page_configuration():
    st.set_page_config(page_title="Pangenome Composition", page_icon=":dna:", layout="wide")
    sys.setrecursionlimit(100000)  # You can adjust the limit as needed

def load_data(matrix_file):
    try:
        if matrix_file is not None:
            data = pd.read_csv(matrix_file, sep='\t', index_col=0)
            return data
    except Exception as e:
        st.error(f"Error loading data: {str(e)}")
    return None

def sample_data(data, threshold):
    try:
        if data is not None:
            if threshold >= data.shape[1]:
                return data
            else:
                np.random.seed(123)
                sampled_columns = np.random.choice(data.columns, threshold, replace=False)
                return data[sampled_columns]
    except Exception as e:
        st.error(f"Error sampling data: {str(e)}")
    return None

def create_pangenome_data(threshold, genes_matrix_df):
    try:
        pangenome_data = pd.DataFrame(columns=["genome", "pangenome", "core", "addition", "unique"])

        # ... (your existing code)

        return pangenome_data
    except Exception as e:
        st.error(f"Error creating pangenome data: {str(e)}")
        return None

def plot_pangenome_analysis(pangenome_data, data_anno, threshold):
    try:
        # ... (your existing code)

        st.pyplot(fig)
        st.pyplot(fig2)
    except Exception as e:
        st.error(f"Error plotting pangenome analysis: {str(e)}")

def pangenome_composition():
    try:
        set_page_configuration()

        # Sidebar for data loading and settings
        st.sidebar.title('Data Loading and Settings')
        matrix_file = st.sidebar.file_uploader('Upload Matrix File (Text File)', type=['txt'])
        sample_button = st.sidebar.button('Sample Data')
        threshold = st.sidebar.slider('Threshold', min_value=1, max_value=500, value=50)

        # Main content
        st.title('Pangenome Analysis')

        if matrix_file is not None:
            data = load_data(matrix_file)
            sampled_data = sample_data(data, threshold)
            
            if sampled_data is not None:
                data_anno = pd.DataFrame({
                    'genes': sampled_data.index,
                    'isolates': sampled_data.sum(axis=1),
                    'group': np.where(
                        (sampled_data.sum(axis=1) < threshold * 0.01) | (sampled_data.sum(axis=1) == 1), 'unique',
                        np.where(sampled_data.sum(axis=1) >= threshold * 0.99, 'core', 'pangenome')
                    )
                })
                data_anno = data_anno[data_anno['isolates'] != 0]

                pangenome_data = create_pangenome_data(threshold, sampled_data)

                # Create a row for the first two plots
                st.write('<style>div.row-widget.stHorizontal {flex-direction: row;}</style>', unsafe_allow_html=True)
                col1, col2 = st.columns(2)
                plot_pangenome_analysis(pangenome_data, data_anno, threshold)
    except Exception as e:
        st.error(f"An unexpected error occurred: {str(e)}")