import os
import sys
import csv
from pprint import pprint as pp

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import Align
from Bio.Align import substitution_matrices
from bayes_opt import BayesianOptimization


def max_score(sequence, aligner):
    """
    Calculate the maximum alignment score of a sequence against itself.
    
    Args:
    - sequence (str): The sequence to evaluate.
    - aligner (Align.PairwiseAligner): The alignment tool.

    Returns:
    - float: Maximum alignment score.
    """
    alignments = aligner.align(sequence, sequence)
    return alignments[0].score


def calculate_similarity(s1, s2, aligner, max_score_s1, max_score_s2):
    """
    Compute the similarity between two sequences.
    
    Args:
    - s1, s2 (str): Sequences to compare.
    - aligner (Align.PairwiseAligner): The alignment tool.
    - max_score_s1, max_score_s2 (float): Max scores for sequences s1 and s2.

    Returns:
    - float: Similarity score.
    """
    alignments = aligner.align(s1, s2)
    actual_score = alignments.score
    return actual_score / ((max_score_s1 + max_score_s2) / 2)


def write_pml_file(pdb, fname, mode):
    """
    Generate a .pml file from a DataFrame based on the provided mode.
    
    Args:
    - pdb (str): PDB name.
    - fname (str): Destination filename.
    - mode (str): Operational mode.
    """
    directory = os.getcwd()
    df = pd.read_csv(f"cosine_similarities_Mode{mode}_updated_ranked_LocalAln_weightedSum.csv")
    df = df.drop_duplicates(subset='sequence_1', keep='first')
    
    with open(fname, "w+") as f:
        f.write(f"load {directory}/{pdb}.pdb\n")
        for sequence in df['sequence_1']:
            f.write(f"select ps. {sequence}\n")
            f.write("color gold, sele\n")


def read_and_calculate_similarity(mode_num, aligner):
    """
    Extract data from CSV and compute sequence similarity.
    
    Args:
    - mode_num (str): The mode number.
    - aligner (Align.PairwiseAligner): The alignment tool.

    Returns:
    - Tuple: Data and its headers.
    """
    data = []
    with open(f"cosine_similarities_Mode{mode_num}.csv", "r") as f:
        reader = csv.reader(f)
        headers = next(reader)
        for row in reader:
            data.append(row)

    headers.append("sequence_similarity")

    # Enhance the data with sequence similarity
    for row in data:
        s1, s2 = row[0], row[1]
        max_score_s1 = max_score(s1, aligner)
        max_score_s2 = max_score(s2, aligner)
        similarity = calculate_similarity(s1, s2, aligner, max_score_s1, max_score_s2)
        row.append(similarity)

    return data, headers


def write_updated_data_to_csv(data, headers, mode_num):
    """
    Save the updated data to a new CSV file.
    
    Args:
    - data (list): Data to save.
    - headers (list): Headers for the data columns.
    - mode_num (str): Mode number.
    """
    with open(f"cosine_similarities_Mode{mode_num}_updated_LocalAln.csv", "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(data)


def calculate_ranks_and_write_to_csv(mode_num):
    """
    Compute ranks and save to CSV.
    
    Args:
    - mode_num (str): Mode number.

    Returns:
    - DataFrame: Data with computed ranks.
    """
    df = pd.read_csv(f"cosine_similarities_Mode{mode_num}_updated_LocalAln.csv")
    df['cosine_rank'] = df['cosine_similarity'].rank(ascending=False)
    df['sequence_rank'] = df['sequence_similarity'].rank(ascending=False)
    df['frobenius_rank'] = df['frobenius_similarity_normed'].rank(ascending=True)
    df['rank_sum'] = df['cosine_rank'] + df['sequence_rank'] + df['frobenius_rank']
    df['rank_sum_pct'] = df['rank_sum'].rank(pct=True)

    df_ranked = df.sort_values('rank_sum_pct')
    df_ranked.to_csv(f"cosine_similarities_Mode{mode_num}_updated_ranked_LocalAln.csv", index=False)

    return df_ranked


def calculate_weighted_sum_whole_and_save(df_ranked, mode_num, w_cosine=0.5, w_frobenius=0.1, w_sw=0.4):
    """
    Compute the weighted sum for sequences and save.
    
    Args:
    - df_ranked (DataFrame): Ranked data.
    - mode_num (str): Mode number.
    - w_cosine (float, optional): Weight for cosine similarity. Default is 0.5.
    - w_frobenius (float, optional): Weight for frobenius similarity. Default is 0.1.
    - w_sw (float, optional): Weight for sequence similarity. Default is 0.4.
    """
    # Normalize weights
    total_weight = w_cosine + w_frobenius + w_sw
    w_cosine /= total_weight
    w_frobenius /= total_weight
    w_sw /= total_weight

    max_length = df_ranked['tile_length'].max()
    dfs = []

    for length in range(6, max_length + 1):
        df_filtered = df_ranked[df_ranked['tile_length'] == length].copy()
        df_filtered['weighted_sum'] = (w_cosine * df_filtered['cosine_similarity'] + 
                                       w_frobenius * df_filtered['frobenius_similarity_normed'] + 
                                       w_sw * df_filtered['sequence_similarity'])
        dfs.append(df_filtered)

    df_combined = pd.concat(dfs)
    df_combined.to_csv(f"cosine_similarities_Mode{mode_num}_updated_ranked_LocalAln_weightedSum_whole.csv", index=False)


def analyze_and_visualize(df):
    """
    Produce statistical descriptions and visualizations.
    
    Args:
    - df (DataFrame): Data to analyze and visualize.
    """
    print("Descriptive Statistics for Sequence Similarity:\n", df['sequence_similarity'].describe())

    # Histogram
    sns.histplot(df['sequence_similarity'], kde=True)
    plt.title('Distribution of Sequence Similarity Scores')
    plt.show()

    # Boxplot
    sns.boxplot(y=df['sequence_similarity'])
    plt.title('Boxplot of Sequence Similarity Scores')
    plt.show()

    # Scatterplot
    sns.scatterplot(data=df, x='cosine_similarity', y='sequence_similarity')
    plt.title('Cosine Similarity vs. Sequence Similarity')
    plt.show()

    # Correlation
    correlation = df['cosine_similarity'].corr(df['sequence_similarity'])
    print(f"\nCorrelation between Cosine Similarity and Sequence Similarity: {correlation:.2f}")

    # Correlation matrix
    corr_matrix = df[['frobenius_similarity_normed', 'cosine_similarity', 'sequence_similarity']].corr()
    print("\nCorrelation Matrix:\n", corr_matrix)

    # Heatmap of correlation matrix
    sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', vmin=-1, vmax=1)
    plt.title('Heatmap of Correlation Matrix')
    plt.show()

    # Threshold analysis suggestion
    if df['sequence_similarity'].median() < 0.50:
        print("\nConsidering the median sequence similarity is below 0.50, you might want to reconsider this threshold.")


def plot_convergence(optimizer):
    """
    Visualize the convergence of the Bayesian optimization process.
    
    Args:
    - optimizer (BayesianOptimization): Optimizer instance.
    """
    y_best = [x['target'] for x in optimizer.res]
    plt.figure(figsize=(10, 5))
    plt.plot(y_best, '-o')
    plt.xlabel('Iterations')
    plt.ylabel('Objective Value')
    plt.title('Convergence Plot')
    plt.grid(True)
    plt.show()


def objective(w_cosine, w_frobenius, w_sw, df):
    """
    Objective function for the Bayesian optimization.
    
    Args:
    - w_cosine, w_frobenius, w_sw (float): Weights.
    - df (DataFrame): Data.

    Returns:
    - float: Objective value.
    """
    total_weight = w_cosine + w_frobenius + w_sw
    w_cosine /= total_weight
    w_frobenius /= total_weight
    w_sw /= total_weight

    df['weighted_sum'] = (w_cosine * df['cosine_similarity'] + 
                          w_frobenius * df['frobenius_similarity_normed'] + 
                          w_sw * df['sequence_similarity'])
    
    # Average of the correlations
    corr_cosine = df['weighted_sum'].corr(df['cosine_similarity'])
    corr_frobenius = df['weighted_sum'].corr(df['frobenius_similarity_normed'])
    corr_sequence = df['weighted_sum'].corr(df['sequence_similarity'])
    avg_corr = (corr_cosine + corr_frobenius + corr_sequence) / 3
    
    return avg_corr


def optimize_weights(df):
    """
    Determine optimal weights using Bayesian Optimization.
    
    Args:
    - df (DataFrame): Data.

    Returns:
    - dict: Optimal normalized weights.
    """
    pbounds = {'w_cosine': (0, 1), 'w_frobenius': (0, 1), 'w_sw': (0, 1)}
    optimizer = BayesianOptimization(
        f=lambda w_cosine, w_frobenius, w_sw: objective(w_cosine, w_frobenius, w_sw, df),
        pbounds=pbounds,
        random_state=1,
    )

    optimizer.maximize(init_points=15, n_iter=100)
    best_params = optimizer.max['params']

    # Normalize weights
    total_weight = sum(best_params.values())
    best_params_normed = {k: v/total_weight for k, v in best_params.items()}

    # Display results
    print("\nOptimal Weights (Normalized):")
    for k, v in best_params_normed.items():
        print(f"{k}: {v:.4f}")

    plot_convergence(optimizer)

    return best_params_normed


def main():
    mode_num = sys.argv[1]

    # Setup aligner
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.mode = 'local'

    # Process sequence similarity
    data, headers = read_and_calculate_similarity(mode_num, aligner)
    write_updated_data_to_csv(data, headers, mode_num)

    # Rank data
    df_ranked = calculate_ranks_and_write_to_csv(mode_num)

    # Weighted sum calculation
    calculate_weighted_sum_whole_and_save(df_ranked, mode_num)
    df_combined = pd.read_csv(f"cosine_similarities_Mode{mode_num}_updated_ranked_LocalAln_weightedSum_whole.csv")

    # Optimize weights
    best_weights = optimize_weights(df_combined)
    pp(best_weights)

    # Recompute weighted sum with optimal weights
    calculate_weighted_sum_whole_and_save(df_ranked, mode_num, **best_weights)
    df_combined = pd.read_csv(f"cosine_similarities_Mode{mode_num}_updated_ranked_LocalAln_weightedSum_whole.csv")

    # Analyze and visualize
    analyze_and_visualize(df_combined)


if __name__ == "__main__":
    main()
