import os
import numpy as np
import pandas as pd
import scanpy as sc
import requests
import tarfile
import matplotlib.pyplot as plt
import configparser
import subprocess
from sklearn.metrics import accuracy_score, f1_score, precision_score
import argparse
import pickle

BASE_DIR      = os.path.dirname(os.path.abspath(__file__))
BASE_DATA_DIR = os.path.join(BASE_DIR, 'dataset')
BASE_IMG_DIR  = os.path.join(BASE_DIR, 'Images')

def parse_add_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument('--dataset_id', type=str, help='Unique identifier of the dataset')

    return parser


def load_anndata(dataset_id):
    dataset_path = os.path.join(BASE_DATA_DIR, dataset_id)
    adata = sc.read_10x_mtx(dataset_path, var_names='gene_symbols', cache=True)
    adata.var_names_make_unique()
    return adata


# def load_anndata(dataset_id):
#     dataset_dir = os.path.join(BASE_DATA_DIR, dataset_id)
    
#     for filename in os.listdir(dataset_dir):
#         if filename.endswith(".tar.gz"):
#             file_path = os.path.join(dataset_dir, filename)
#             with tarfile.open(file_path, 'r:gz') as tar:
#                 tar.extractall(dataset_dir)
#                 print(f"Extracted {filename} in {dataset_dir}")


def preprocessing(dataset_id,adata,Base_IMG_DIR):
    # Step 2: Preprocessing
    # Show genes with the highest fraction of counts in each single cell
    #root_dir = os.path.dirname(os.path.abspath(__file__))
    #images_dir = os.path.join(root_dir, 'Images')
    #os.makedirs(images_dir, exist_ok=True)
    highest_expr_genes_path = os.path.join(Base_IMG_DIR, f'{dataset_id}_highest_expr_genes.png')
    # Plot the figure
    sc.pl.highest_expr_genes(adata, n_top=20, show=False)
    plt.title('Highest Expressed genes')
    
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
    violin_plot_path = os.path.join(Base_IMG_DIR, f'{dataset_id}_violin_plot.png')
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
    plt.savefig(violin_plot_path, dpi=300)
    plt.close()   
    
    # Remove cells with too many mitochondrial genes or too many total counts
    scatter_plot_path = os.path.join(Base_IMG_DIR, f'{dataset_id}_scatter_plot.png')
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
    highly_variable_gene_path = os.path.join(Base_IMG_DIR, f'{dataset_id}_highly_variable_gene.png')
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3.0, min_disp=0.5)
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



def run_pca(adata,dataset_id,images_dir):
    sc.tl.pca(adata, svd_solver='arpack', n_comps=60)
    sc.pl.pca(adata, color='CST3')
    PCA_path = os.path.join(images_dir, f'{dataset_id}_PCA.png')
    plt.savefig(PCA_path, dpi=300)
    plt.close()
    print('data before annotation\n',adata.obs)
    return adata



cluster_to_cell_type_dict = {}
def single_cell_analysis(adata,dataset_id,num_neighbors_values, num_pcs_values, images_dir):
    '*****************COMPUTING AND EMBEDDING NEIGHBOURHOOD GRAPH***********'
    pca_neighbors_dict = {}
    for i, num_pcs in enumerate(num_pcs_values):
        for j, num_neighbors in enumerate(num_neighbors_values):
            adata_copy = adata.copy()
            sc.pp.neighbors(adata_copy, n_neighbors=num_neighbors, n_pcs=num_pcs)
            sc.tl.leiden(adata_copy)
            sc.tl.paga(adata_copy)
            sc.pl.paga(adata_copy, plot=False)
            sc.tl.umap(adata_copy, init_pos='paga')
            sc.tl.umap(adata_copy)
            title_str = f"PCA-{num_pcs}_Neighbors-{num_neighbors}"
            # Plot UMAP with the custom title
            sc.pl.umap(adata_copy, color=['CST3', 'NKG7'], use_raw=False, title=title_str)
            
            # Save the UMAP figure with the custom title
            umap_path = os.path.join(images_dir, f'{dataset_id}_UMAP_PCA-{num_pcs}_Neighbors-{num_neighbors}.png')
            plt.title(f'UMAP for PCA-{num_pcs}_Neighbors-{num_neighbors} in dataset_{dataset_id}')
            plt.savefig(umap_path, dpi=300)
            plt.close() 
            
            '*****************CLustering NEIGHBOURHOOD GRAPH***********'
            
            sc.tl.leiden(adata_copy, key_added='leiden')
            sc.pl.umap(adata_copy, color=['leiden', 'CST3'])
            UMAP_path_leiden = os.path.join(images_dir, f'{dataset_id}_UMAP_WITH_LEIDEN_PCA-{num_pcs}_Neighbors-{num_neighbors}.png')
            plt.title(f'Umap with Leiden CLustering for PCA-{num_pcs}_Neighbors-{num_neighbors} in dataset_{dataset_id}')
            plt.savefig(UMAP_path_leiden, dpi=300)
            plt.close()
            
            '*****************Finding Marker Genes ***********'
            sc.tl.rank_genes_groups(adata_copy, 'leiden', method='t-test')
            sc.pl.rank_genes_groups(adata_copy, n_genes=25, sharey=False)
            marker_genes_path = os.path.join(images_dir, f'{dataset_id}_Finding_marker_genes_PCA-{num_pcs}_Neighbors-{num_neighbors}.png')
            plt.title(f'Marker genes for PCA-{num_pcs}_Neighbors-{num_neighbors} in dataset_{dataset_id}')
            plt.savefig(marker_genes_path, dpi=300)
            plt.close()
            
            result = adata_copy.uns["rank_genes_groups"]
            groups = result["names"].dtype.names
            
            '*****************creating SCSA result***********'
            
            df = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'logfoldchanges','scores','pvals']})
            df.to_csv(f"{dataset_id}_diff_exp_result_PCA-{num_pcs}_Neighbors-{num_neighbors}.csv")
            
            # Generate input and output file names dynamically based on PCA and number of neighbors
            input_file_name = f"{dataset_id}_diff_exp_result_PCA-{num_pcs}_Neighbors-{num_neighbors}.csv"
            output_file_name = f"{dataset_id}_scsa_result_PCA-{num_pcs}_Neighbors-{num_neighbors}.txt"
            
            # Run the script using subprocess with dynamically generated input and output file names
            script_path = "./run_scsa.sh"
            try:
                subprocess.run(["bash", script_path, "-i", input_file_name, "-o", output_file_name], check=True)
                print("run_scsa.sh executed successfully.")
            except subprocess.CalledProcessError as e:
                print("Error executing run_scsa.sh:", e)
            
            '****************** Annotation*********************'
            df = pd.read_csv(output_file_name, sep='\t')
            # Group the data by "Cluster" and find the cell type with the highest Z-score in each group
            highest_zscores = df.groupby('Cluster')['Z-score'].idxmax()
            print ('**********************************************')
            highest_zscores.fillna('Unknown', inplace=True)
            print("Index of highest zscore after replacement:", highest_zscores.index)
            # Extract the corresponding cell types for the highest Z-scores
            df.reset_index(inplace=True, drop=True)
            cell_types_with_highest_zscores = df.loc[highest_zscores.index, 'Cell Type'].tolist()

            #cell_types_with_highest_zscores = df.loc[highest_zscores, 'Cell Type'].tolist()
            print ('cell type with highest Z-Score ' , cell_types_with_highest_zscores)
            
            #  Create a dictionary of cluster value (number) and corresponding identified cell type. 
            #  We will add a column to adata with cell types for each cluster number.
            cluster_to_cell_type = {cluster_num: cell_type for cluster_num, cell_type in enumerate(cell_types_with_highest_zscores)}
            print('Cluster_to_cell_type ', cluster_to_cell_type)
            cluster_to_cell_type_dict[(num_pcs, num_neighbors)] = cluster_to_cell_type

            print ('printing cluster_to_cell_type_dict',cluster_to_cell_type_dict)
            print("***********************")
            print('copied adata', adata_copy)
            print("***********************")
        
            adata_copy.obs['leiden'] = adata_copy.obs['leiden'].astype(int)
            # Map cluster numbers to cell type 
            adata_copy.obs['cell_types'] = adata_copy.obs['leiden'].map(cluster_to_cell_type)
            #unique cell types 
            print('unique cell types', adata_copy.obs['cell_types'].unique())
            sc.pl.umap(adata_copy, color=["leiden", "cell_types"])
            final_cell_type_data = os.path.join(images_dir, f'{dataset_id}_final_cell_type_PCA-{num_pcs}_Neighbors-{num_neighbors}.png')
            plt.title(f'SCSA cell type annotation for PCA-{num_pcs}_Neighbors-{num_neighbors}for_{dataset_id}')
            plt.savefig(final_cell_type_data, dpi=300)
            plt.close()
            '******************Calculating and plotting Evaluation matrics*********************'
            if num_pcs == 40 and num_neighbors == 10:
            # Extract the golden standard only when PCA is 40 and number of neighbors is 10
                golden_standard = adata_copy.obs['cell_types']
                if np.nan in golden_standard.unique():
                    golden_standard = golden_standard.astype('str')
                    golden_standard = golden_standard.replace('nan', 'Unknown')
            
            predicted_cell_types = adata_copy.obs['cell_types']
            if np.nan in predicted_cell_types.unique():
                predicted_cell_types = predicted_cell_types.astype('str')
                predicted_cell_types = predicted_cell_types.replace('nan', 'Unknown')
            
            # Store the PCA and number of neighbors values along with the corresponding cell types
            pca_neighbors_dict[(num_pcs, num_neighbors)] = predicted_cell_types
            

    print("Golden Standard data types:", golden_standard.dtype)
    print("Predicted Cell Types data types:", predicted_cell_types.dtype)
    print('copied adata\n',adata_copy.obs) 
    # Initialize lists to store evaluation metrics for each combination
    accuracies = []
    precisions = []
    f1_scores = []

    # Exclude the point where PCA=40 and num_neighbors=10
    excluded_combination = (40, 10)
    # Initialize a separate list to store num_pcs_values for evaluation metrics
    eval_num_pcs_values = []
    eval_num_neighbour_values = []
    for (num_pcs, num_neighbors), predicted_cell_types in pca_neighbors_dict.items():
        if (num_pcs, num_neighbors) != excluded_combination:
            # Calculate accuracy by comparing predicted cell types with the golden standard
            print("="*100)
            print(type(golden_standard), type(predicted_cell_types))
            print(golden_standard)
            print(predicted_cell_types)
            print("="*100)
            accuracy = accuracy_score(golden_standard, predicted_cell_types)
            precision = precision_score(golden_standard, predicted_cell_types, average='weighted')
            f1 = f1_score(golden_standard, predicted_cell_types, average='weighted')

            # Print or store the accuracy, precision, and F1 score for the current combination
            print(f"Accuracy for PCA={num_pcs} and num_neighbors={num_neighbors}: {accuracy}")
            print(f"Precision for PCA={num_pcs} and num_neighbors={num_neighbors}: {precision}")
            print(f"F1 score for PCA={num_pcs} and num_neighbors={num_neighbors}: {f1}")
            
            # Append the metrics to the respective lists
            accuracies.append(accuracy)
            precisions.append(precision)
            f1_scores.append(f1)
            eval_num_pcs_values.append(num_pcs)
            eval_num_neighbour_values.append(num_neighbors)
    # Plotting the evaluation metrics against the number of principal components
    plt.figure(figsize=(10, 6))
    plt.plot(eval_num_neighbour_values, accuracies, marker='o', label='Accuracy')
    plt.plot(eval_num_neighbour_values, precisions, marker='o', label='Precision')
    plt.plot(eval_num_neighbour_values, f1_scores, marker='o', label='F1 Score')

    # Add labels and title
    plt.xlabel('Number of Neighbours')
    plt.ylabel('Score')
    plt.title(f'Change of Evaluation Metrics with Number of Neighbours for {dataset_id}')
    plt.legend()
    plt.grid(True)

    # Save the plot
    accuracy_plot_path = os.path.join(images_dir, f'{dataset_id}_evaluation_metrics_vs_NumOfNeighbours.png')
    plt.savefig(accuracy_plot_path, dpi=300)
    plt.close()
    # Plotting the evaluation metrics against the number of principal components
    plt.figure(figsize=(10, 6))
    plt.plot(eval_num_pcs_values, accuracies, marker='o', label='Accuracy')
    plt.plot(eval_num_pcs_values, precisions, marker='o', label='Precision')
    plt.plot(eval_num_pcs_values, f1_scores, marker='o', label='F1 Score')

    # Add labels and title
    plt.xlabel('Number of Principal Components (PCA)')
    plt.ylabel('Score')
    plt.title(f'Change of Evaluation Metrics with PCA for {dataset_id}')
    plt.legend()
    plt.grid(True)

    # Save the plot
    accuracy_plot_path = os.path.join(images_dir, f'{dataset_id}_evaluation_metrics_vs_pca.png')
    plt.savefig(accuracy_plot_path, dpi=300)
    plt.close()
    return adata_copy


def save_adata(adata_copy, file_path):
    """
    Save the AnnData object to a file.

    Parameters:
        adata (AnnData): The AnnData object to save.
        file_path (str): The file path to save the AnnData object to.
    """
    if adata_copy is not None:
        adata_copy.write(file_path)
        print(f"AnnData object saved to {file_path}")
    else:
        print(f"AnnData object saved to {file_path}")




def read_config(config_file_name):
    config = configparser.ConfigParser()
    config.read(config_file_name)
    return config

def main():
    print(BASE_DIR)
    print(BASE_IMG_DIR)
    parser = argparse.ArgumentParser(description='Pipeline Script')
    parser = parse_add_args(parser)
    args   = parser.parse_args()
    
    dataset_id = args.dataset_id
    adata = load_anndata(dataset_id)
    adata = run_pca(adata, dataset_id, images_dir=BASE_IMG_DIR)
    # Step 2: Preprocessing
    adata = preprocessing(dataset_id,adata,Base_IMG_DIR=BASE_IMG_DIR)
    #single cell analysis
    config = read_config('config2.ini') 
    # Read num_neighbors_values and num_pcs_values from config
    num_neighbors_values = [int(x) for x in config['Parameters']['num_neighbors_values'].split(',')]
    num_pcs_values = [int(x) for x in config['Parameters']['num_pcs_values'].split(',')]

    adata = single_cell_analysis(adata, dataset_id, num_neighbors_values, num_pcs_values, images_dir=BASE_IMG_DIR)
    
    # # step 5: saving results before annotation
    results_path = f"{dataset_id}_processed.h5ad"
    save_adata(adata, results_path)

if __name__ == "__main__":
    main()
