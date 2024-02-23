import os
import numpy as np
import pandas as pd
import scanpy as sc
import requests
import tarfile
import matplotlib.pyplot as plt


def download_data(url, data_fname, save_dir,extracted_file_path):
    os.makedirs(save_dir, exist_ok=True)
    file_path = os.path.join(save_dir, data_fname)
    
    res = requests.get(url)
    with open(file_path, 'wb') as file:
        file.write(res.content)
    
    # Extract the downloaded file
    with tarfile.open(file_path, 'r:gz') as tar:
        tar.extractall(save_dir)
    
    
    
    print(f"Data downloaded and extracted to {extracted_file_path}")
    return extracted_file_path


def preprocessing(extracted_file_path, min_mean, max_mean, images_dir):
    # Step 1: Read in the count matrix into an AnnData object
    adata = sc.read_10x_mtx(
        extracted_file_path,
        var_names='gene_symbols',
        cache=True
    )
    adata.var_names_make_unique()
    
    
    
    # Step 2: Preprocessing
    # Show genes with the highest fraction of counts in each single cell
    #root_dir = os.path.dirname(os.path.abspath(__file__))
    #images_dir = os.path.join(root_dir, 'Images')
    #os.makedirs(images_dir, exist_ok=True)
    highest_expr_genes_path = os.path.join(images_dir, 'highest_expr_genes.png')
    # Plot the figure
    sc.pl.highest_expr_genes(adata, n_top=20, show=False)
    
    # Save the figure
    plt.savefig(highest_expr_genes_path, dpi=300)
    plt.close()    
    

    
    
    # Basic filtering
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Assemble information about mitochondrial genes for quality control
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Violin plot of computed quality measures
    violin_plot_path = os.path.join(images_dir, 'violin_plot.png')
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
    plt.savefig(violin_plot_path, dpi=300)
    plt.close()   
    
    
    
    # Remove cells with too many mitochondrial genes or too many total counts
    scatter_plot_path = os.path.join(images_dir, 'scatter_plot.png')
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
    plt.savefig(scatter_plot_path, dpi=300)
    plt.close()  
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    
    # Total-count normalize the data matrix to 10,000 reads per cell
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    # Logarithmize the data
    sc.pp.log1p(adata)
    
    # Identify highly-variable genes
    highly_variable_gene_path = os.path.join(images_dir, 'highly_variable_gene.png')
    sc.pp.highly_variable_genes(adata, min_mean=min_mean, max_mean=max_mean, min_disp=0.5)
    sc.pl.highly_variable_genes(adata)
    plt.savefig(highly_variable_gene_path, dpi=300)
    plt.close()
    
    
    # Set the .raw attribute of the AnnData object to the normalized and logarithmized raw gene expression
    adata.raw = adata
    
    # Actually do the filtering
    adata = adata[:, adata.var.highly_variable]
    
    # Regress out effects of total counts per cell and the percentage of mitochondrial genes
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    
    # Scale the data to unit variance
    sc.pp.scale(adata, max_value=10)
    
    
    # Return the processed AnnData object
    return adata




def run_pca(adata, n_pcs,images_dir):
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_pcs)
    sc.pl.pca(adata, color='CST3')
    PCA_path = os.path.join(images_dir, 'PCA.png')
    plt.savefig(PCA_path, dpi=300)
    plt.close()
    return adata
    
    
def computing_and_embedding_neighbourhood_graphs(adata, num_neighbors_values, num_pcs_values, images_dir):
    num_subplots = len(num_neighbors_values) * len(num_pcs_values)
    fig, axes = plt.subplots(len(num_pcs_values), len(num_neighbors_values), figsize=(15, 5 * len(num_pcs_values)))

    for i, num_pcs in enumerate(num_pcs_values):
        for j, num_neighbors in enumerate(num_neighbors_values):
            sc.pp.neighbors(adata, n_neighbors=num_neighbors, n_pcs=num_pcs)
            sc.tl.leiden(adata)
            sc.tl.paga(adata)
            sc.pl.paga(adata, plot=False)
            sc.tl.umap(adata, init_pos='paga')
            sc.tl.umap(adata)
            title_str = f"PCA-{num_pcs}_Neighbors-{num_neighbors}"
            sc.pl.umap(adata, color=['CST3'], use_raw=False, title=title_str, ax=axes[i, j])
            
     #Save the UMAP figure with the custom title
    umap_path = os.path.join(images_dir,'UMAP.png')
    plt.savefig(umap_path, dpi=300)
    plt.close()
    return adata
    
# def ClusteringNeighbourhoodGraphs(adata,images_dir):
#     sc.tl.leiden(adata)
#     sc.pl.umap(adata, color=['leiden', 'CST3', 'NKG7'])
#     UMAP_path_leiden = os.path.join(images_dir, 'UMAP_WITH_LEIDEN.png')
#     plt.savefig(UMAP_path_leiden, dpi=300)
#     plt.close()
#     return adata

'''
methods to try: 
t-test
wilcoxon (Wilcoxon Rank Sum Test)
t-test_overestim_var (t-test with Welch's Correction)
ANOVA (Analysis of Variance)
kruskal (Kruskal-Wallis Test)
fisher (Fisher's Exact Test)
'''

# def finding_marker_genes(num_genes,method, adata,images_dir): 
#     sc.tl.rank_genes_groups(adata, 'leiden', method=method)
#     sc.pl.rank_genes_groups(adata, n_genes=num_genes, sharey=False)
#     marker_genes_path = os.path.join(images_dir, 'Finding_marker_genes.png')
#     plt.savefig(marker_genes_path, dpi=300)
#     plt.close()
#     save_adata(adata, "processed_adata.h5ad")
    
    
def save_adata(adata, file_path):
    """
    Save the AnnData object to a file.

    Parameters:
        adata (AnnData): The AnnData object to save.
        file_path (str): The file path to save the AnnData object to.
    """
    if adata is not None:
        adata.write(file_path)
        print(f"AnnData object saved to {file_path}")
    else:
        print(f"AnnData object saved to {file_path}")

##################################################################################################################################################################################


import configparser

def read_config(file_path):
    config = configparser.ConfigParser()
    config.read(file_path)
    return config

def main():
    # Step 1: Data Downloading
    url = "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
    data_fname = "pbmc3k_filtered_gene_bc_matrices.tar.gz"
    extracted_file_path = 'data/filtered_gene_bc_matrices/hg19/'
    
    root_dir = os.path.dirname(os.path.abspath(__file__))
   
    images_dir = os.path.join(root_dir, 'Images')
    save_dir = os.path.join(root_dir, 'data')
    os.makedirs(images_dir, exist_ok=True)
    extracted_file_path = download_data(url, data_fname, save_dir,extracted_file_path)

    # Step 2: Preprocessing
    config = read_config('config.ini')
    min_mean = float(config['Preprocessing']['min_mean'])
    max_mean = float(config['Preprocessing']['max_mean'])
    adata = preprocessing(extracted_file_path, min_mean, max_mean, images_dir)

    # # Step 3: PCA
    n_pcs = int(config['PCA']['n_pcs'])
    adata = run_pca(adata, n_pcs, images_dir)

    # # Step 4: Computing and embedding neighbourhood graphs
    
    num_neighbors_values = [5, 10, 12]  
    num_pcs_values = [20, 30]  
    adata = computing_and_embedding_neighbourhood_graphs(adata, num_neighbors_values, num_pcs_values, images_dir)
    
    #num_pcs = int(config['NeighbourhoodGraphs']['num_pcs'])
    #adata = computing_and_embedding_neighbourhood_graphs(adata, num_neighbors_values, num_pcs, images_dir)

    # # # Step 5: Clustering neighbourhood graph
    # adata = ClusteringNeighbourhoodGraphs(adata, images_dir)

    # # # Step 6: finding marker genes
    # num_genes = int(config['FindingMarkerGenes']['num_genes'])
    # method = config['FindingMarkerGenes']['method']
    # adata = finding_marker_genes(num_genes, method, adata, images_dir)

    results_path = "pbmc3k_processed.h5ad"
    save_adata(adata, results_path)
    
    
if __name__ == "__main__":
    main()

