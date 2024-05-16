import os
import pickle
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(BASE_DIR, 'result')
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

    # Rearrange the DataFrame to have all accuracies, precisions, and F1 scores together
    summary_table_acc = summary_table.filter(like='Accuracy')
    summary_table_prec = summary_table.filter(like='Precision')
    summary_table_f1 = summary_table.filter(like='F1 Score')
    rearranged_summary_table = pd.concat([summary_table_acc, summary_table_prec, summary_table_f1], axis=1)

    # Plot the summary table
    plt.figure(figsize=(16, 10))
    sns.heatmap(rearranged_summary_table, annot=False, fmt=".3f", cmap="YlGnBu")
    plt.title('Summary Table', fontsize=20)
    plt.xlabel('Metrics', fontsize=16)
    plt.ylabel('Datasets', fontsize=16)
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'summary_table.png'))
    plt.show()

if __name__ == "__main__":
    generate_summary_table()
