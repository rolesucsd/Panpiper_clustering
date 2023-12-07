import streamlit as st
from ete3 import Tree, TreeStyle, NodeStyle
import pandas as pd
import seaborn as sns
from io import BytesIO
import base64
import tempfile
import threading

def load_data(metadata_file, tree_file):
    metadata = pd.read_csv(metadata_file, index_col=0, sep=None, engine='python')
    tree = Tree(tree_file.read().decode("utf-8"))
    return metadata, tree

def render_tree_image(metadata, tree, metadata_category):
    with threading.Lock():
        midpoint_node = tree.get_midpoint_outgroup()

        # Reroot the tree at the midpoint node
        tree.set_outgroup(midpoint_node)

        ts = TreeStyle()
        ts.mode = "c"
        ts.show_leaf_name = False

        # Assuming metadata is your DataFrame and metadata_category is the column of interest
        unique_values = metadata[metadata_category].unique()

        # Use seaborn to generate a color palette based on the number of unique values
        num_unique_values = len(unique_values)
        color_palette = sns.color_palette("husl", n_colors=num_unique_values)

        # Create a dictionary mapping unique values to colors using the #FFFFFF format
        metadata_colors = {value: f"#{int(color[0]*255):02X}{int(color[1]*255):02X}{int(color[2]*255):02X}" for value, color in zip(unique_values, color_palette)}

        for n in tree.traverse():
            nstyle = NodeStyle()
            # Hide node circles
            if n.is_leaf():
                n_metadata = metadata.at[n.name, metadata_category]
                node_color = metadata_colors.get(n_metadata, "#000000")  # Default to black if no color is specified
                nstyle["fgcolor"] = node_color
                nstyle["size"] = 15
                n.set_style(nstyle)

        # Create a temporary file to save the tree
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as temp_file:
            temp_file_name = temp_file.name
            tree.render(temp_file_name, tree_style=ts)

        temp_file.close()

        return temp_file_name

def phylogeny():
    # Streamlit app title and description
    st.title("Phylogenetic Tree Visualization with ETE Toolkit")
    st.write("This app allows you to visualize a circular phylogenetic tree with colored tips based on metadata.")

    # Upload files
    metadata_file = st.sidebar.file_uploader("Upload metadata TXT file", type=["txt"])
    tree_file = st.sidebar.file_uploader("Upload Newick tree file", type=["nwk", "newick"])

    # Sidebar for user input
    st.sidebar.title("User Input")

    if metadata_file and tree_file:
        metadata, tree = load_data(metadata_file, tree_file)
        metadata_category = st.sidebar.selectbox("Select metadata category for coloring tips", metadata.columns.tolist())
        image_path = render_tree_image(metadata, tree, metadata_category)

        # Display the tree image using Streamlit
        st.image(image_path, use_column_width=True)
    else:
        st.warning("Please upload both a metadata TXT file and a Newick tree file.")
