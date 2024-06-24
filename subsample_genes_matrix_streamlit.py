import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

st.set_page_config(page_title="Pangenome Composition", page_icon=":dna:", layout="wide")
sys.setrecursionlimit(100000)  # You can adjust the limit as needed

# Function to load data from file
def load_data(matrix_file):
    if matrix_file is not None:
        data = pd.read_csv(matrix_file, sep='\t', index_col=0)
        return data
    return None

# Function to sample data based on threshold
def sample_data(data, threshold):
    if data is not None:
        if threshold >= data.shape[1]:
            return data
        else:
            np.random.seed(123)
            sampled_columns = np.random.choice(data.columns, threshold, replace=False)
            return data[sampled_columns]
    return None

import pandas as pd

# Function to create pangenome data
def create_pangenome_data(threshold, genes_matrix_df):
    pangenome_data = pd.DataFrame(columns=["genome", "pangenome", "core", "addition", "unique"])
    data_list = []

    size = genes_matrix_df.shape[1]
    d1 = genes_matrix_df.iloc[:, 0]

    initial_data = {
        "genome": 1, 
        "pangenome": len(d1[d1 != 0]), 
        "core": len(d1[d1 == 1]), 
        "addition": len(d1[d1 != 0]), 
        "unique": len(d1[d1 != 0])
    }
    data_list.append(initial_data)

    for i in range(2, threshold + 1):
        genome = i
        d1 = d1 + genes_matrix_df.iloc[:, i - 1].apply(pd.to_numeric)
        pangenome = len(d1[d1 != 0])
        coregenome = len(d1[d1 >= i * 0.99])
        addition = pangenome - data_list[-1]["pangenome"]
        unique = len(d1[d1 == 1])

        data = {
            "genome": genome, 
            "pangenome": pangenome, 
            "core": coregenome, 
            "addition": addition, 
            "unique": unique
        }
        data_list.append(data)

    pangenome_data = pd.DataFrame(data_list)

    return pangenome_data


# Main function
def pangenome_composition():
    # Sidebar for data loading and settings
    st.set_option('deprecation.showfileUploaderEncoding', False)
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
            # Pangenome analysis
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

        with col1: 
            # Pangenome Size Analysis Plot
            st.subheader('Pangenome Size Analysis')
            fig, ax = plt.subplots()
            sns.lineplot(data=pangenome_data, x='genome', y='core', label='Core Genes')
            sns.lineplot(data=pangenome_data, x='genome', y='pangenome', label='Pangenome')
            sns.lineplot(data=pangenome_data, x='genome', y='unique', label='Unique Genes')
            st.pyplot(fig)

            # Addition Plot
            x_ticks = range(0, threshold, 10)

            st.subheader('Addition Plot')
            fig, ax = plt.subplots()
            sns.histplot(x='isolates', hue='group', data=data_anno, bins=threshold)
            ax.set_xticks(x_ticks)

            st.pyplot(fig)

        with col2: 
            # Calculate the threshold for row sum
            threshold = 0.01 * len(sampled_data.columns)
            # Filter rows based on the threshold
            subsampled_data = sampled_data[sampled_data.sum(axis=1) > threshold]
            st.subheader('Subsampled Heatmap')
            heatmap = sns.clustermap(subsampled_data, cmap="YlGnBu", xticklabels=False, yticklabels=False, annot=False, figsize=(10, 20), cbar=False, method='average', metric='euclidean')
            st.pyplot(heatmap)
    else:
        st.warning('Please upload a matrix file and sample the data.')
