import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import pdist, squareform
from skbio.stats.distance import permanova
from scipy.cluster.hierarchy import linkage, fcluster, cut_tree
from skbio import DistanceMatrix
import random
from sklearn.cluster import DBSCAN
from sklearn.manifold import MDS
from sklearn.manifold import TSNE
import os  # Add this line to import the os module
import umap
import matplotlib.pyplot as plt
import streamlit as st
from io import BytesIO
import base64

@st.cache  # Use Streamlit caching for data loading
def load_data(mash_file, group_file):
    mash = pd.read_csv(mash_file, sep="\t", index_col=0)
    mash.index = [os.path.splitext(os.path.basename(row))[0] for row in mash.index]
    mash.columns = [os.path.splitext(os.path.basename(col))[0] for col in mash.columns]
    group = pd.read_csv(group_file, sep="\t", encoding='latin-1')
    return mash, group

def perform_clustering(mash, clustering_height):
    dist = mash
    square_dist = squareform(dist)
    linkage_matrix = linkage(square_dist, method='ward')

    max_d = clustering_height * np.max(linkage_matrix[:, 2])
    clusters = cut_tree(linkage_matrix, height=max_d)

    cluster_labels = [str(i + 1) for i in clusters.flatten()]
    # Plot clusters
    cluster_labels = [str(label).strip('[]') for label in clusters]  # Remove "[" and "]" characters
    phylogroup = pd.DataFrame({'Phylogroup': cluster_labels, 'Sample': mash.index})

    return cluster_labels, phylogroup

def perform_permanova(mash, group, selection):
    group_filtered = group[group[selection].isin([0, 1])]
    common_ids = set(mash.index).intersection(group_filtered['Sample'])
    mash_filtered = mash.loc[common_ids, common_ids]
    group_filtered = group_filtered[group_filtered['Sample'].isin(common_ids)]

    sample_order = mash_filtered.index.tolist()
    grouping = group_filtered.set_index('Sample').loc[sample_order, selection].tolist()

    mash_df = mash_filtered.apply(pd.to_numeric, errors='coerce')
    mash_df = np.ascontiguousarray(mash_df.values)
    mash_dm = DistanceMatrix(mash_df, ids=common_ids)

    result = permanova(distance_matrix=mash_dm, grouping=grouping, permutations=10000)
    return result

def plot_manual_clustering(group, selection, pca_df, tsne_df, name1, name2):
    unique_selection = group[selection].unique()
    random.shuffle(unique_selection)
    selection_colors = {group: '#' + ''.join(random.choices('0123456789ABCDEF', k=6)) for group in unique_selection}

    pca_df = pd.merge(pca_df, group, left_on='SampleName', right_on='Sample', how='inner')
    fig, ax = plt.subplots()
    sns.scatterplot(x=name1, y=name2, hue=selection, palette=selection_colors, data=pca_df, ax=ax)
    ax.legend(title='Cluster', loc='center left', bbox_to_anchor=(1, 0.5))  # Adjust the legend location

    st.pyplot(fig)  # Display the figure in Streamlit

    tsne_df = pd.merge(tsne_df, group, left_on='SampleName', right_on='Sample', how='inner')
    fig, ax = plt.subplots()
    sns.scatterplot(x='TSNE1', y='TSNE2', hue=selection, palette=selection_colors, data=tsne_df, ax=ax)
    ax.legend(title='Cluster', loc='center left', bbox_to_anchor=(1, 0.5))  # Adjust the legend location
    st.pyplot(fig)  # Display the figure in Streamlit


def pca_and_tsne(mash, cluster_labels):
    pca = PCA(n_components=2)
    pca.fit(mash)
    pca_scores = pca.transform(mash)
    name1 = f'PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)'
    name2 = f'PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)'
    pca_df = pd.DataFrame(pca_scores, columns=[name1, name2])
    pca_df['Cluster'] = cluster_labels
    pca_df['SampleName'] = mash.index  # Assuming index contains sample names

    tsne = TSNE(n_components=2)
    embedding = tsne.fit_transform(mash)
    tsne_df = pd.DataFrame(embedding, columns=['TSNE1', 'TSNE2'])
    tsne_df['Cluster'] = cluster_labels
    tsne_df['SampleName'] = mash.index  # Assuming index contains sample names

    return pca_df, tsne_df, name1, name2


def plot_pca_and_tsne(pca_df, tsne_df, phylogroup, name1, name2):
    unique_phylogroups = phylogroup['Phylogroup'].unique()
    random.shuffle(unique_phylogroups)
    phylogroup_colors = {group: '#' + ''.join(random.choices('0123456789ABCDEF', k=6)) for group in unique_phylogroups}

    fig, ax = plt.subplots()
    sns.scatterplot(x=name1, y=name2, hue='Cluster', palette=phylogroup_colors, data=pca_df, ax=ax)
    ax.legend(title='Cluster', loc='center left', bbox_to_anchor=(1, 0.5))  # Adjust the legend location
    st.pyplot(fig)  # Display the figure in Streamlit

    fig, ax = plt.subplots()
    sns.scatterplot(x='TSNE1', y='TSNE2', hue='Cluster', palette=phylogroup_colors, data=tsne_df, ax=ax)
    ax.legend(title='Cluster', loc='center left', bbox_to_anchor=(1, 0.5))  # Adjust the legend location
    st.pyplot(fig)  # Display the figure in Streamlit
    return 0

def download_link(object_to_download, download_filename, download_link_text):
    if isinstance(object_to_download, pd.DataFrame):
        object_to_download = object_to_download.to_csv(index=False)

    # Some strings <-> bytes conversions necessary here
    b64 = base64.b64encode(object_to_download.encode()).decode()
    return f'<a href="data:file/txt;base64,{b64}" download="{download_filename}">{download_link_text}</a>'


def main():
    # Set page title and font
    st.set_page_config(page_title='Phylogroup Selection and Grouping', page_icon=':dna:')
    st.title('Phylogroup Selection and Grouping')
    st.markdown('<style>h1{font-family: "MonoLisa", monospace; font-size: 32px !important;}</style>', unsafe_allow_html=True)

    # Create a sidebar for user input
    st.sidebar.title('User Input')
    mash_file = st.sidebar.file_uploader('Upload MASH File (TSV)', type=['tsv'])
    group_file = st.sidebar.file_uploader('Upload Group Data (Text File)', type=['txt'])
    clustering_height = st.sidebar.slider('Clustering Height', 0.1, 1.0, 0.44, 0.01)
    selection = st.sidebar.text_input('Manual Clustering Selection (e.g., "Isolation Source")')

    phylogroup = None  # Define phylogroup

    if mash_file is not None and group_file is not None:
        st.write('Data loaded successfully.')

        mash, group = load_data(mash_file, group_file)
        cluster_labels, phylogroup = perform_clustering(mash, clustering_height)
        pca_df, tsne_df, name1, name2 = pca_and_tsne(mash, cluster_labels)
        plot_pca_and_tsne(pca_df, tsne_df, phylogroup, name1, name2)

        result = perform_permanova(mash, group, selection)
        st.header('PerMANOVA Result')
        st.write(result)
        plot_manual_clustering(group, selection, pca_df, tsne_df, name1, name2)

    if st.button("Download Phylogroup Data"):
        if phylogroup is not None:
            # Provide a download link for the phylogroup data
            tmp_download_link = download_link(phylogroup, 'phylogroup.txt', 'Click here to download phylogroup.txt')
            st.markdown(tmp_download_link, unsafe_allow_html=True)

if __name__ == '__main__':
    main()