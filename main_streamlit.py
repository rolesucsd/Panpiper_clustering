import streamlit as st
from subsample_genes_matrix_streamlit import pangenome_composition
from phylogroup_cluster_streamlit import phylogroup
from diff_abund_streamlit import abundance

def main():
    st.sidebar.title("Navigation")
    app_selector = st.sidebar.radio("Go to", ["Pangenome Composition", "Phylogrouping", "Differential Prevalence"])

    if app_selector == "Pangenome Composition":
        pangenome_composition()

    elif app_selector == "Phylogrouping":
        phylogroup()

    elif app_selector == "Differential Prevalence":
        abundance()

if __name__ == "__main__":
    main()
