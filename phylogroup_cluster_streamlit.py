import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from scipy.spatial.distance import pdist, squareform
from skbio.stats.distance import permanova
from scipy.cluster.hierarchy import linkage, cut_tree
from skbio import DistanceMatrix
import matplotlib.pyplot as plt
import streamlit as st
from io import BytesIO
import base64
import random

def load_data(mash_file, group_file):
    try:
        mash = pd.read_csv(mash_file, sep="\t", index_col=0, header=0)
        group = pd.read_csv(group_file, sep=None, engine='python')
        return mash, group
    except Exception as e:
        st.error(f"Error loading data: {str(e)}")
        return None, None

def extract_sample_id(index):
    parts = index.split('/')
    return parts[-2] if len(parts) >= 2 else index

def perform_clustering(mash, clustering_height):
    try:
        square_dist = squareform(mash)
        linkage_matrix = linkage(square_dist, method='ward')

        max_d = clustering_height * np.max(linkage_matrix[:, 2])
        clusters = cut_tree(linkage_matrix, height=max_d)

        cluster_labels = [str(i + 1) for i in clusters.flatten()]
        phylogroup = pd.DataFrame({'Phylogroup': cluster_labels, 'Sample': mash.index})

        return cluster_labels, phylogroup
    except Exception as e:
        st.error(f"Error performing clustering: {str(e)}")
        return None, None

def perform_permanova(mash, group, selection):
    try:
        common_ids = list(set(mash.index).intersection(group['Sample']))
        mash_filtered = mash.loc[common_ids, common_ids]
        group_filtered = group[group['Sample'].isin(common_ids)].drop_duplicates(subset="Sample", keep="first")

        if mash_filtered.empty or group_filtered.empty:
            st.error("Filtered data is empty. Check your input data.")
            return None

        sample_order = mash_filtered.index.tolist()
        grouping = group_filtered.set_index('Sample').loc[sample_order, selection].tolist()

        mash_df = mash_filtered.apply(pd.to_numeric, errors='coerce')
        mash_df = np.ascontiguousarray(mash_df.values)
        mash_dm = DistanceMatrix(mash_df, ids=common_ids)

        if len(grouping) != len(mash_dm.ids):
            st.error("Number of samples in grouping and distance matrix do not match.")
            return None

        result = permanova(distance_matrix=mash_dm, grouping=grouping, permutations=10000)
        return result
    except Exception as e:
        st.error(f"Error performing PerMANOVA: {str(e)}")
        return None

def plot_manual_clustering(group, selection, pca_df, tsne_df, name1, name2):
    try:
        unique_selection = group[selection].unique()
        random.shuffle(unique_selection)
        selection_colors = {group: '#' + ''.join(random.choices('0123456789ABCDEF', k=6)) for group in unique_selection}

        pca_df = pd.merge(pca_df, group, left_on='SampleName', right_on='Sample', how='inner')
        fig, ax = plt.subplots()
        sns.scatterplot(x=name1, y=name2, hue=selection, palette=selection_colors, data=pca_df, ax=ax)
        ax.legend(title='Cluster', loc='center left', bbox_to_anchor=(1, 0.5))

        st.pyplot(fig)

        tsne_df = pd.merge(tsne_df, group, left_on='SampleName', right_on='Sample', how='inner')
        fig, ax = plt.subplots()
        sns.scatterplot(x='TSNE1', y='TSNE2', hue=selection, palette=selection_colors, data=tsne_df, ax=ax)
        ax.legend(title='Cluster', loc='center left', bbox_to_anchor=(1, 0.5))

        st.pyplot(fig)
    except Exception as e:
        st.error(f"Error plotting manual clustering: {str(e)}")

def pca_and_tsne(mash, cluster_labels):
    try:
        pca = PCA(n_components=2)
        pca.fit(mash)
        pca_scores = pca.transform(mash)
        name1 = f'PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)'
        name2 = f'PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)'
        pca_df = pd.DataFrame(pca_scores, columns=[name1, name2])
        pca_df['Cluster'] = cluster_labels
        pca_df['SampleName'] = mash.index

        tsne = TSNE(n_components=2)
        embedding = tsne.fit_transform(mash)
        tsne_df = pd.DataFrame(embedding, columns=['TSNE1', 'TSNE2'])
        tsne_df['Cluster'] = cluster_labels
        tsne_df['SampleName'] = mash.index

        return pca_df, tsne_df, name1, name2
    except Exception as e:
        st.error(f"Error performing PCA and t-SNE: {str(e)}")
        return None, None, None, None

def plot_pca_and_tsne(pca_df, tsne_df, phylogroup, name1, name2):
    try:
        unique_phylogroups = phylogroup['Phylogroup'].unique()
        random.shuffle(unique_phylogroups)
        phylogroup_colors = {group: '#' + ''.join(random.choices('0123456789ABCDEF', k=6)) for group in unique_phylogroups}

        fig, ax = plt.subplots()
        sns.scatterplot(x=name1, y=name2, hue='Cluster', palette=phylogroup_colors, data=pca_df, ax=ax)
        ax.legend(title='Cluster', loc='center left', bbox_to_anchor=(1, 0.5))

        st.pyplot(fig)

        fig, ax = plt.subplots()
        sns.scatterplot(x='TSNE1', y='TSNE2', hue='Cluster', palette=phylogroup_colors, data=tsne_df, ax=ax)
        ax.legend(title='Cluster', loc='center left', bbox_to_anchor=(1, 0.5))

        st.pyplot(fig)
    except Exception as e:
        st.error(f"Error plotting PCA and t-SNE: {str(e)}")

def download_link(object_to_download, download_filename, download_link_text):
    try:
        if isinstance(object_to_download, pd.DataFrame):
            object_to_download = object_to_download.to_csv(index=False)

        b64 = base64.b64encode(object_to_download.encode()).decode()
        st.markdown(f'<a href="data:file/txt;base64,{b64}" download="{download_filename}">{download_link_text}</a>', unsafe_allow_html=True)
    except Exception as e:
        st.error(f"Error creating download link: {str(e)}")

def phylogroup():
    try:
        st.title('Phylogroup Selection and Grouping')
        st.markdown('<style>h1{font-family: "MonoLisa", monospace; font-size: 32px !important;}</style>', unsafe_allow_html=True)

        st.sidebar.title('User Input')
        mash_file = st.sidebar.file_uploader('Upload MASH File (TSV)', type=['tsv'])
        group_file = st.sidebar.file_uploader('Upload Group Data (Text File)', type=['txt'])
        clustering_height = st.sidebar.slider('Clustering Height', 0.1, 1.0, 0.44, 0.01)

        phylogroup = None

        if mash_file and group_file:
            mash, group = load_data(mash_file, group_file)
            group_columns = list(group.columns)
        else:
            st.warning('Please upload a mash distance matrix.')
            return

        selection = st.sidebar.selectbox("Clustering Selection", group_columns)
        group = group.rename(columns={group.columns[0]: "Sample"})
        mash.index = mash.index.map(extract_sample_id)
        mash.columns = mash.columns.map(extract_sample_id)

        cluster_labels, phylogroup = perform_clustering(mash, clustering_height)
        if cluster_labels and phylogroup is not None:
            pca_df, tsne_df, name1, name2 = pca_and_tsne(mash, cluster_labels)
            plot_pca_and_tsne(pca_df, tsne_df, phylogroup, name1, name2)

            result = perform_permanova(mash, group, selection)
            if result is not None:
                st.header('PerMANOVA Result')
                st.write(result)
                plot_manual_clustering(group, selection, pca_df, tsne_df, name1, name2)

            if st.button("Download Phylogroup Data"):
                download_link(phylogroup, 'phylogroup.txt', 'Click here to download phylogroup.txt')

    except Exception as e:
        st.error(f"An unexpected error occurred: {str(e)}")

if __name__ == '__main__':
    phylogroup()
