import os
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from pandas.plotting import table
import numpy as np

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(BASE_DIR, 'result')
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, 'dataset')

def generate_summary_table():
    summary_table = pd.DataFrame()
    dataset_dirs = [d for d in os.listdir(DATA_DIR) if os.path.isdir(os.path.join(DATA_DIR, d))]

    for dataset_id in dataset_dirs:
        # Load evaluation results for the current dataset
        evaluation_results_file = os.path.join(RESULTS_DIR, f'{dataset_id}_evaluation_results.pkl')
        with open(evaluation_results_file, 'rb') as f:
            evaluation_results = pickle.load(f)

        # Extract evaluation metrics
        accuracies = evaluation_results['accuracy']
        precisions = evaluation_results['precision']
        f1_scores = evaluation_results['f1_score']

        # Populate summary table with evaluation metrics
        for (num_pcs, num_neighbors), accuracy in accuracies.items():
            col_prefix = f"PCA{num_pcs}_NB{num_neighbors}"
            summary_table.at[dataset_id, f'Accuracy_{col_prefix}'] = accuracy
            summary_table.at[dataset_id, f'Precision_{col_prefix}'] = precisions[(num_pcs, num_neighbors)]
            summary_table.at[dataset_id, f'F1 Score_{col_prefix}'] = f1_scores[(num_pcs, num_neighbors)]

    # Save the table image
    table_image_path = os.path.join(RESULTS_DIR, 'summary_table.png')
    fig, ax = plt.subplots(figsize=(16, 10))
    ax.axis('off')
    table(ax, summary_table, loc='center', cellLoc='center', colWidths=[0.15] * len(summary_table.columns),
          cellColours=plt.cm.Greens(np.full(summary_table.shape, 0.5)))
    plt.savefig(table_image_path, bbox_inches='tight', pad_inches=0.05)
    plt.show()
if __name__ == "__main__":
    generate_summary_table()
