# Refactoring and Improving Imports

# Standard Libraries
import os
import sys
import csv
from pprint import pprint as pp
from ast import literal_eval
from itertools import product

# Third-party Libraries
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.PDB import PDBParser
from scipy.spatial.distance import cosine
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm



# Utility Functions

def read_fasta(file_path):
    """
    Reads sequences from a FASTA file.
    
    Parameters:
    - file_path (str): Path to the FASTA file.
    
    Returns:
    - list of str: List of sequences from the FASTA file.
    """
    with open(file_path, "r") as fasta_file:
        sequences = [str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]
    return sequences


def cosine_similarity(A, B):
    """
    Calculates the cosine similarity between two arrays.
    
    Parameters:
    - A (array-like): First array.
    - B (array-like): Second array.
    
    Returns:
    - float: Cosine similarity between A and B.
    """
    dot_product = np.dot(A, B)
    norm_A = np.linalg.norm(A)
    norm_B = np.linalg.norm(B)
    return dot_product / (norm_A * norm_B)


def frobenius_similarity(A, B):
    """
    Calculates the Frobenius norm of the difference between two matrices.
    
    Parameters:
    - A (matrix): First matrix.
    - B (matrix): Second matrix.
    
    Returns:
    - float: Frobenius norm of the difference between A and B.
    """
    diff = A - B
    return np.linalg.norm(diff, ord='fro')


def normalize_data(data):
    """
    Performs min-max normalization on the data.
    
    Parameters:
    - data (array-like): Data to be normalized.
    
    Returns:
    - array-like: Normalized data.
    """
    min_val = np.min(data)
    max_val = np.max(data)
    return (data - min_val) / (max_val - min_val)


def pearson_correlation(A, B):
    """
    Calculates the Pearson correlation coefficient between two arrays.
    
    Parameters:
    - A (array-like): First array.
    - B (array-like): Second array.
    
    Returns:
    - tuple: Pearson correlation coefficient and p-value.
    """
    return pearsonr(A, B)


def get_tiles(ANM_crossCorrelation_mat, sequence, lower_tile_length):
    """
    Extracts submatrices (tiles) and corresponding sequence tiles from a given matrix.
    
    Parameters:
    - ANM_crossCorrelation_mat: Numpy matrix from which tiles are extracted.
    - sequence: Corresponding sequence data.
    - lower_tile_length: Minimum length of tiles to extract.

    Returns:
    - List of dictionaries containing tiles and metadata.
    """
    N = len(ANM_crossCorrelation_mat)
    tiles = [
        {
            'tile': ANM_crossCorrelation_mat[start:end, start:end],
            'sequence': sequence[start:end],
            'length': L,
            'start_index': start + 1,
            'end_index': end,
            'center_index': np.median(np.arange(start, end)) + 1
        }
        for L in range(N, lower_tile_length-1, -1)
        for start in range(N - L + 1)
        for end in [start + L]
    ]
    return tiles


def get_residue_range(pdb_file, chain_id):
    """
    Returns the range of residue numbers for a given chain in a PDB file.
    
    Parameters:
    - pdb_file: Path to the PDB file.
    - chain_id: ID of the chain for which the residue range is to be found.

    Returns:
    - Tuple of (minimum residue number, maximum residue number).
    """
    # Initialize a PDBParser object
    parser = PDBParser()

    # Read the structure
    structure = parser.get_structure('pdb_structure', pdb_file)

    # Get the specific chain
    chain = structure[0][chain_id]

    # Extract residue numbers and determine min and max
    residue_numbers = [residue.get_id()[1] for residue in chain]
    return min(residue_numbers), max(residue_numbers)


def group_tiles_by_length(tiles):
    """
    Groups tiles by their length.
    
    Parameters:
    - tiles: List of tiles.

    Returns:
    - Dictionary with keys as tile lengths and values as lists of tiles of that length.
    """
    tiles_by_length = {}
    for index, tile in tiles.iterrows():
        tiles_by_length.setdefault(tile['length'], []).append(tile)
    
    return tiles_by_length


def write_similarity_to_csv(writer, tile_1, tile_2, length, similarity, frobenius_sim, within_protein=False):
    """
    Writes the similarity data to a CSV file. Can handle both within-protein and between-proteins comparisons.
    
    Parameters:
    - writer: CSV writer object.
    - tile_1: Tile from the first protein or first tile from the same protein.
    - tile_2: Tile from the second protein or second tile from the same protein.
    - length: Length of the tile.
    - similarity: Calculated cosine similarity.
    - frobenius_sim: Calculated Frobenius similarity.
    - within_protein: Boolean flag to indicate if the comparison is within the same protein.
    """
    if within_protein:
        writer.writerow([
            tile_1['sequence'], tile_1['start_index'], tile_1['end_index'], tile_1['center_index'],
            tile_2['sequence'], tile_2['start_index'], tile_2['end_index'], tile_2['center_index'],
            length, similarity, frobenius_sim
        ])
    else:
        writer.writerow([
            tile_1['sequence'], tile_2['sequence'], 
            tile_1['start_index'], tile_1['end_index'], tile_1['center_index'], 
            tile_2['start_index'], tile_2['end_index'], tile_2['center_index'], 
            length, similarity, frobenius_sim
        ])


def calculate_cosine_similarities_between_proteins(tiles_protein_1, tiles_protein_2, modeNum):
    """
    Calculates the cosine similarities between tiles of two proteins and saves the results to a CSV file.
    
    Parameters:
    - tiles_protein_1: List of tiles from the first protein.
    - tiles_protein_2: List of tiles from the second protein.
    - modeNum: Mode number (used for naming the output CSV file).
    """
    # Group tiles by their length
    tiles_by_length_1 = group_tiles_by_length(tiles_protein_1)
    tiles_by_length_2 = group_tiles_by_length(tiles_protein_2)

    with open(f'cosine_similarities_Mode{modeNum}.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["sequence_1", "sequence_2", "start_index_1", "end_index_1", "center_index_1", 
                         "start_index_2", "end_index_2", "center_index_2", "tile_length", 
                         "cosine_similarity", "frobenius_similarity"])

        # Iterate over lengths present in both proteins
        for length in set(tiles_by_length_1.keys()) & set(tiles_by_length_2.keys()):
            for tile_2 in tqdm(tiles_by_length_2[length], desc=f"Processing length {length}"):
                for tile_1 in tiles_by_length_1[length]:
                    flat_tile_1 = tile_1['tile'].flatten()
                    flat_tile_2 = tile_2['tile'].flatten()
                    similarity = cosine_similarity(flat_tile_1, flat_tile_2)
                    frobenius_sim = frobenius_similarity(tile_1['tile'], tile_2['tile'])
                    write_similarity_to_csv(writer, tile_1, tile_2, length, similarity, frobenius_sim)

    print("Finished comparing all matching tile lengths.")


def calculate_cosine_similarities_within_protein(tiles_protein, modeNum):
    """
    Calculates the cosine similarities between tiles within the same protein and saves the results to a CSV file.
    
    Parameters:
    - tiles_protein: List of tiles from the protein.
    - modeNum: Mode number (used for naming the output CSV file).
    """
    # Group tiles by their length
    tiles_by_length = group_tiles_by_length(tiles_protein)

    with open(f'cosine_similarities_Mode{modeNum}_within_protein.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["sequence_1", "start_index_1", "end_index_1", "center_index_1", 
                         "sequence_2", "start_index_2", "end_index_2", "center_index_2", 
                         "tile_length", "cosine_similarity", "frobenius_similarity"])

        # Iterate over all lengths
        for length, tiles in tiles_by_length.items():
            for tile_1 in tqdm(tiles, desc=f"Processing length {length}"):
                for tile_2 in tiles:
                    flat_tile_1 = tile_1['tile'].flatten()
                    flat_tile_2 = tile_2['tile'].flatten()
                    similarity = cosine_similarity(flat_tile_1, flat_tile_2)
                    frobenius_sim = frobenius_similarity(tile_1['tile'], tile_2['tile'])
                    write_similarity_to_csv(writer, tile_1, tile_2, length, similarity, frobenius_sim, within_protein=True)

    print("Finished comparing all matching tile lengths within the protein.")



# Here's a refactored version that breaks down the main function:

def read_and_process_cross_correlation_matrix(csv_dir, pdb_id, anm_mode_number, first_res, last_res):
    """
    Reads and processes the cross-correlation matrix.
    
    Parameters:
    - csv_dir: Directory containing the CSV file.
    - pdb_id: ID of the protein.
    - anm_mode_number: Mode number for ANM.
    - first_res: First residue index.
    - last_res: Last residue index.

    Returns:
    - Processed cross-correlation matrix.
    """
    csv_file = f"{pdb_id}_CrossCorr_Mode{anm_mode_number}.csv"
    csv_path = os.path.join(csv_dir, csv_file)
    
    ANM_crossCorrelation_mat = pd.read_csv(csv_path, header=None).values
    return ANM_crossCorrelation_mat[first_res:last_res + 1, first_res:last_res + 1]




# def process_protein_tiles(ANM_crossCorrelation_mat, sequence, pdb_id):
#     if not os.path.exists(f'{pdb_id}_tiles.pkl'):
#         tiles_protein = get_tiles(ANM_crossCorrelation_mat, sequence, 6)
#         prot_tile_df = pd.DataFrame.from_dict(tiles_protein)
#         prot_tile_df = prot_tile_df.sort_values(by=['length', 'center_index']).reset_index(drop=True)
#         # prot_tile_df.to_pickle(f'{pdb_id}_tiles.pkl')
#     else:
#         prot_tile_df = pd.read_pickle(f'{pdb_id}_tiles.pkl')
    
#     return prot_tile_df


def process_protein_tiles(ANM_crossCorrelation_mat, sequence, pdb_id):
    """
    Processes the protein tiles. If the pickle file exists, it reads it. 
    Otherwise, it computes the tiles and saves them as a pickle file.
    
    Parameters:
    - ANM_crossCorrelation_mat: The cross-correlation matrix.
    - sequence: Protein sequence.
    - pdb_id: ID of the protein.

    Returns:
    - DataFrame of protein tiles.
    """
    pickle_file = f'{pdb_id}_tiles.pkl'
    
    if not os.path.exists(pickle_file):
        tiles_protein = get_tiles(ANM_crossCorrelation_mat, sequence, 6)
        prot_tile_df = pd.DataFrame.from_dict(tiles_protein).sort_values(by=['length', 'center_index']).reset_index(drop=True)
    else:
        prot_tile_df = pd.read_pickle(pickle_file)
    
    return prot_tile_df




# def process_output_csv(anm_mode_number, protein_length):
#     tile_df = pd.read_csv(f'cosine_similarities_Mode{anm_mode_number}.csv')
#     tile_df['frobenius_similarity_normed'] = normalize_data(tile_df['frobenius_similarity'])
#     tile_df.to_csv(f'cosine_similarities_Mode{anm_mode_number}.csv', index=False)

#     threshold = 1e-9
#     condition = (tile_df['tile_length'] >= 6) & (tile_df['tile_length'] <= protein_length) & (np.abs(tile_df['cosine_similarity'] - 1) > threshold)
#     filtered_tile_df = tile_df[condition]
    
#     grouped = filtered_tile_df.groupby(['center_index_1', 'tile_length'])
#     average_cosine_similarity = grouped['cosine_similarity'].mean()
#     counts = grouped.size()

#     result = pd.DataFrame({'average_cosine_similarity': average_cosine_similarity, 'counts': counts}).reset_index()
#     result.rename(
#         {
#             'center_index_1': 'Tile_Center_Index',
#             'tile_length': 'Tile_Length',
#             'average_cosine_similarity': 'Average_Cosine_Similarity',
#             'counts': 'Count'
#         },
#         axis=1,
#         inplace=True
#     )
#     return result


def process_output_csv(anm_mode_number, protein_length):
    """
    Processes the output CSV file. Calculates and adds the normalized frobenius similarity.
    Filters the data based on tile length and cosine similarity conditions.

    Parameters:
    - anm_mode_number: Mode number for ANM.
    - protein_length: Length of the protein sequence.

    Returns:
    - DataFrame containing the filtered average cosine similarities.
    """
    csv_file = f'cosine_similarities_Mode{anm_mode_number}.csv'
    tile_df = pd.read_csv(csv_file)
    tile_df['frobenius_similarity_normed'] = normalize_data(tile_df['frobenius_similarity'])
    tile_df.to_csv(csv_file, index=False)

    threshold = 1e-9
    condition = (tile_df['tile_length'] >= 6) & (tile_df['tile_length'] <= protein_length) & (np.abs(tile_df['cosine_similarity'] - 1) > threshold)
    filtered_tile_df = tile_df[condition]
    
    grouped = filtered_tile_df.groupby(['center_index_1', 'tile_length'])
    average_cosine_similarity = grouped['cosine_similarity'].mean()
    counts = grouped.size()

    result = pd.DataFrame({'average_cosine_similarity': average_cosine_similarity, 'counts': counts}).reset_index()
    result.rename(columns={
        'center_index_1': 'Tile_Center_Index',
        'tile_length': 'Tile_Length',
        'average_cosine_similarity': 'Average_Cosine_Similarity',
        'counts': 'Count'
    }, inplace=True)
    
    return result








# The main function is then simplified as follows:

def main():
    args = sys.argv
    csv_dir = os.path.join("../", "ANM_modes")
    print(args)
    
    if len(args) == 4:
        pdb_id, anm_mode_number, chain_id = args[1:4]

        pdb_file = f"{pdb_id}.chain{chain_id}.nowa.pdb"
        first_res, last_res = get_residue_range(pdb_file, chain_id)
        
        ANM_crossCorrelation_mat = read_and_process_cross_correlation_matrix(csv_dir, pdb_id, anm_mode_number, first_res, last_res)
        concat_sequence = read_fasta(os.path.join(os.getcwd(), f"{pdb_id}.nowa.fasta"))
        prot_tile_df = process_protein_tiles(ANM_crossCorrelation_mat, concat_sequence[0], pdb_id)
        
        protein_length = len(concat_sequence[0])
        calculate_cosine_similarities_within_protein(prot_tile_df, anm_mode_number)
        
        result = process_output_csv(anm_mode_number, protein_length)
        result.to_csv(f'filtered_average_cosine_similarities_Mode{anm_mode_number}_within_protein.csv', index=False)
    
    elif len(args) == 5 or len(args) == 6:
        pdb_id_one, pdb_id_two, anm_mode_number, chain_id_one, chain_id_two = args[1:6]

        # Process for the first protein
        pdb_file_one = f"{pdb_id_one}.chain{chain_id_one}.nowa.pdb"
        first_res_one, last_res_one = get_residue_range(pdb_file_one, chain_id_one)
        ANM_crossCorrelation_mat_1 = read_and_process_cross_correlation_matrix(csv_dir, pdb_id_one, anm_mode_number, first_res_one, last_res_one)
        concat_sequence_A = read_fasta(os.path.join(os.getcwd(), f"{pdb_id_one}.chain{chain_id_one}.nowa.fasta"))
        prot_tile_df_one = process_protein_tiles(ANM_crossCorrelation_mat_1, concat_sequence_A[0], pdb_id_one)
        
        # Process for the second protein
        pdb_file_two = f"{pdb_id_two}.nowa.pdb"
        first_res_two, last_res_two = get_residue_range(pdb_file_two, chain_id_two)
        ANM_crossCorrelation_mat_2 = read_and_process_cross_correlation_matrix(csv_dir, pdb_id_two, anm_mode_number, first_res_two, last_res_two)
        concat_sequence_B = read_fasta(os.path.join(os.getcwd(), f"{pdb_id_two}.nowa.fasta"))
        prot_tile_df_two = process_protein_tiles(ANM_crossCorrelation_mat_2, concat_sequence_B[0], pdb_id_two)
        
        # Calculate the cosine similarities between the two proteins
        calculate_cosine_similarities_between_proteins(prot_tile_df_one, prot_tile_df_two, anm_mode_number)
        
        min_protein_length = min(len(concat_sequence_A[0]), len(concat_sequence_B[0]))
        result = process_output_csv(anm_mode_number, min_protein_length)
        result.to_csv(f'filtered_average_cosine_similarities_Mode{anm_mode_number}_between_chain_{chain_id_one}_and_{chain_id_two}.csv', index=False)


if __name__ == "__main__":
    main()







































# # Refactoring the next set of utility functions and logic

# def get_tiles(ANM_crossCorrelation_mat, sequence, lower_tile_length):
#     """
#     Extracts tiles (submatrices) from the cross-correlation matrix.
    
#     Parameters:
#     - ANM_crossCorrelation_mat (matrix): The cross-correlation matrix.
#     - sequence (str): Protein sequence.
#     - lower_tile_length (int): Minimum length of the tile.
    
#     Returns:
#     - list of dict: List of tiles with details.
#     """
#     N = len(ANM_crossCorrelation_mat)
#     tiles = []
#     for L in range(N, lower_tile_length-1, -1):
#         for start_index in range(N - L + 1):
#             end_index = start_index + L
#             sub_mat = ANM_crossCorrelation_mat[start_index:end_index, start_index:end_index]
#             sequence_tile = sequence[start_index:end_index]
#             tile_center_index = np.median(np.arange(start_index, end_index, 1))
#             tiles.append({
#                 'tile': sub_mat,
#                 'sequence': sequence_tile,
#                 'length': L,
#                 'start_index': start_index + 1,
#                 'end_index': end_index,
#                 'center_index': tile_center_index + 1
#             })
#     return tiles

# def get_residue_range(pdb_file, chain_id):
#     """
#     Retrieves the minimum and maximum residue numbers for a given chain ID in a PDB file.
    
#     Parameters:
#     - pdb_file (str): Path to the PDB file.
#     - chain_id (str): Chain ID to extract residue range.
    
#     Returns:
#     - tuple: Minimum and maximum residue numbers.
#     """
#     parser = PDBParser()
#     structure = parser.get_structure('pdb_structure', pdb_file)
#     chain = structure[0][chain_id]
#     min_res = min(residue.get_id()[1] for residue in chain)
#     max_res = max(residue.get_id()[1] for residue in chain)
#     return min_res, max_res


# # Refactoring functions related to calculating similarities

# def group_tiles_by_length(tiles_df):
#     """
#     Groups tiles by their length.
    
#     Parameters:
#     - tiles_df (DataFrame): DataFrame containing tile information.
    
#     Returns:
#     - dict: Tiles grouped by their length.
#     """
#     tiles_by_length = {}
#     for _, tile in tiles_df.iterrows():
#         tiles_by_length.setdefault(tile['length'], []).append(tile)
#     return tiles_by_length

# def write_to_csv(filename, headers, rows):
#     """
#     Writes data to a CSV file.
    
#     Parameters:
#     - filename (str): Name of the CSV file.
#     - headers (list of str): Headers for the CSV file.
#     - rows (list of list): Data rows to be written to the CSV.
    
#     Returns:
#     - None
#     """
#     with open(filename, 'w', newline='') as file:
#         writer = csv.writer(file)
#         writer.writerow(headers)
#         writer.writerows(rows)

# def calculate_similarities_and_write_to_csv(tiles_1, tiles_2, filename, modeNum):
#     """
#     Calculate cosine and frobenius similarities between tiles and write to CSV.
    
#     Parameters:
#     - tiles_1 (dict): Tiles from the first protein.
#     - tiles_2 (dict): Tiles from the second protein.
#     - filename (str): Name of the CSV file.
#     - modeNum (int): Mode number.
    
#     Returns:
#     - None
#     """
#     # Define CSV headers based on whether it's within or between proteins
#     if tiles_1 is tiles_2:
#         headers = ["sequence_1", "start_index_1", "end_index_1", "center_index_1", 
#                    "sequence_2", "start_index_2", "end_index_2", "center_index_2", 
#                    "tile_length", "cosine_similarity", "frobenius_similarity"]
#     else:
#         headers = ["sequence_1", "sequence_2", "start_index_1", "end_index_1", "center_index_1", 
#                    "start_index_2", "end_index_2", "center_index_2", 
#                    "tile_length", "cosine_similarity", "frobenius_similarity"]
    
#     rows = []
#     # Iterate over lengths present in both sets of tiles
#     for length in set(tiles_1.keys()) & set(tiles_2.keys()):
#         for i in tqdm(range(len(tiles_2[length])), desc=f"Processing length {length}"):
#             for j in range(len(tiles_1[length])):
#                 tile_1 = tiles_1[length][j]
#                 tile_2 = tiles_2[length][i]
#                 flat_tile_1 = tile_1['tile'].flatten()
#                 flat_tile_2 = tile_2['tile'].flatten()
#                 similarity = cosine_similarity(flat_tile_1, flat_tile_2)
#                 frobenius_sim = frobenius_similarity(tile_1['tile'], tile_2['tile'])
#                 row_data = [tile_1['sequence'], tile_1['start_index'], tile_1['end_index'], tile_1['center_index'],
#                             tile_2['sequence'], tile_2['start_index'], tile_2['end_index'], tile_2['center_index'],
#                             length, similarity, frobenius_sim]
#                 if tiles_1 is not tiles_2:
#                     row_data.insert(1, tile_2['sequence'])
#                 rows.append(row_data)
    
#     # Write data to CSV
#     write_to_csv(filename, headers, rows)

# def calculate_cosine_similarities_between_proteins(tiles_protein_1, tiles_protein_2, modeNum):
#     tiles_by_length_1 = group_tiles_by_length(tiles_protein_1)
#     tiles_by_length_2 = group_tiles_by_length(tiles_protein_2)
#     filename = f'cosine_similarities_Mode{modeNum}.csv'
#     calculate_similarities_and_write_to_csv(tiles_by_length_1, tiles_by_length_2, filename, modeNum)
#     print("Finished comparing all matching tile lengths.")

# def calculate_cosine_similarities_within_protein(tiles_protein, modeNum):
#     tiles_by_length = group_tiles_by_length(tiles_protein)
#     filename = f'cosine_similarities_Mode{modeNum}_within_protein.csv'
#     calculate_similarities_and_write_to_csv(tiles_by_length, tiles_by_length, filename, modeNum)
#     print("Finished comparing all matching tile lengths within the protein.")
    
    
# def prepare_tiles_for_protein(csv_path, fasta_path, min_tile_length=6):
#     """
#     Prepares tiles for a given protein based on its cross-correlation matrix and sequence.
    
#     Parameters:
#     - csv_path (str): Path to the CSV file containing the cross-correlation matrix.
#     - fasta_path (str): Path to the FASTA file containing the protein sequence.
#     - min_tile_length (int): Minimum length of the tile (default is 6).
    
#     Returns:
#     - DataFrame: DataFrame containing tile information.
#     """
#     # Read in the cross-correlation matrix
#     ANM_crossCorrelation_mat = pd.read_csv(csv_path, header=None).values
    
#     # Read in the protein sequence
#     concat_sequence = read_fasta(fasta_path)
    
#     # Get tiles from the protein
#     tiles_protein = get_tiles(ANM_crossCorrelation_mat, concat_sequence[0], min_tile_length)
#     prot_tile_df = pd.DataFrame.from_dict(tiles_protein)
#     prot_tile_df = prot_tile_df.sort_values(by=['length', 'center_index'])
#     prot_tile_df = prot_tile_df.reset_index(drop=True)
    
#     return prot_tile_df


# def process_within_protein_data(prot_tile_df, mode_number):
#     """
#     Calculate cosine similarities within the protein, format and save the output data.
    
#     Parameters:
#     - prot_tile_df (DataFrame): DataFrame containing tile information.
#     - mode_number (str): ANM mode number.
    
#     Returns:
#     - DataFrame: DataFrame containing the processed tile data.
#     """
#     # Calculate the Cosine Similarities
#     calculate_cosine_similarities_within_protein(prot_tile_df, mode_number)
    
#     # Reading in and formatting the outputted CSV from the previous function
#     tile_df = pd.read_csv(f'cosine_similarities_Mode{mode_number}_within_protein.csv')
#     tile_df['frobenius_similarity_normed'] = normalize_data(tile_df['frobenius_similarity']) 
#     tile_df.to_csv(f'cosine_similarities_Mode{mode_number}_within_protein.csv', index=False)
    
#     threshold = 1e-9  # Define a small threshold. Adjust this as needed.
#     protein_length = len(prot_tile_df['sequence'][0])
#     condition = (tile_df['tile_length'] >= 6) & (tile_df['tile_length'] <= protein_length / 2) & (np.abs(tile_df['cosine_similarity'] - 1) > threshold)
#     filtered_tile_df = tile_df[condition]
    
#     # Calculate the average cosine similarity and counts for each tile length
#     grouped = filtered_tile_df.groupby(['center_index_1', 'tile_length'])
#     average_cosine_similarity = grouped['cosine_similarity'].mean()
#     counts = grouped.size()

#     # Combine the average_cosine_similarity and counts into a single DataFrame
#     result = pd.DataFrame({'average_cosine_similarity': average_cosine_similarity, 'counts': counts}).reset_index()
#     result.rename(
#         {
#             'center_index_1': 'Tile_Center_Index',
#             'tile_length': 'Tile_Length',
#             'average_cosine_similarity': 'Average_Cosine_Similarity',
#             'counts': 'Count'
#         },
#         axis=1,
#         inplace=True)
    
#     return result


# def process_between_proteins_data(prot_tile_df_one, prot_tile_df_two, mode_number, chain_id_one, chain_id_two):
#     """
#     Calculate cosine similarities between two proteins, format and save the output data.
    
#     Parameters:
#     - prot_tile_df_one (DataFrame): DataFrame containing tile information for the first protein.
#     - prot_tile_df_two (DataFrame): DataFrame containing tile information for the second protein.
#     - mode_number (str): ANM mode number.
#     - chain_id_one (str): Chain ID for the first protein.
#     - chain_id_two (str): Chain ID for the second protein.
    
#     Returns:
#     - DataFrame: DataFrame containing the processed tile data.
#     """
#     # Calculate the Cosine Similarities
#     calculate_cosine_similarities_between_proteins(prot_tile_df_one, prot_tile_df_two, mode_number)
    
#     # Reading in and formatting the outputted CSV from the previous function
#     tile_df = pd.read_csv(f'cosine_similarities_Mode{mode_number}.csv')
#     tile_df['frobenius_similarity_normed'] = normalize_data(tile_df['frobenius_similarity']) 
#     tile_df.to_csv(f'cosine_similarities_Mode{mode_number}.csv', index=False)
    
#     threshold = 1e-9  # Define a small threshold. Adjust this as needed.
#     protein_length_A = len(prot_tile_df_one['sequence'][0])
#     protein_length_B = len(prot_tile_df_two['sequence'][0])
#     condition = (tile_df['tile_length'] >= 6) & (tile_df['tile_length'] <= min(protein_length_A, protein_length_B)) & (np.abs(tile_df['cosine_similarity'] - 1) > threshold)
#     filtered_tile_df = tile_df[condition]
    
#     # Calculate the average cosine similarity and counts for each tile length
#     grouped = filtered_tile_df.groupby(['center_index_1', 'tile_length'])
#     average_cosine_similarity = grouped['cosine_similarity'].mean()
#     counts = grouped.size()

#     # Combine the average_cosine_similarity and counts into a single DataFrame
#     result = pd.DataFrame({'average_cosine_similarity': average_cosine_similarity, 'counts': counts}).reset_index()
#     result.rename(
#         {
#             'center_index_1': 'Tile_Center_Index',
#             'tile_length': 'Tile_Length',
#             'average_cosine_similarity': 'Average_Cosine_Similarity',
#             'counts': 'Count'
#         },
#         axis=1,
#         inplace=True)
    
#     return result



        
# def main():
#     args = sys.argv
#     csv_dir = os.path.join("../", "ANM_modes")
    
#     if len(args) == 4:
#         pdb_id = args[1].upper()
#         anm_mode_number = args[2]
#         chain_id = args[3]
        
#         csv_file = f"{pdb_id}_CrossCorr_Mode{anm_mode_number}.csv"
#         csv_path = os.path.join(csv_dir, csv_file)
#         fasta_path = os.path.join(os.getcwd(), f"{args[1]}.nowa.fasta")
        
#         if not os.path.exists(f'{args[1]}_tiles.pkl'):
#             prot_tile_df = prepare_tiles_for_protein(csv_path, fasta_path)
#             # Save the dataframe to a CSV file
#             prot_tile_df.to_csv(f'{args[1]}_tiles.csv', index=False)
#             # Save the dataframe to a pickle file
#             prot_tile_df.to_pickle(f'{args[1]}_tiles.pkl')
#         else:
#             # If the file does exist, read the pickle file into prot_tile_df
#             prot_tile_df = pd.read_pickle(f'{args[1]}_tiles.pkl')
        
#         # Process data within the protein
#         result = process_within_protein_data(prot_tile_df, anm_mode_number)
#         # Output the result to CSV
#         result.to_csv(f'filtered_average_cosine_similarities_Mode{anm_mode_number}_within_protein.csv', index=False)
        

#     elif len(args) == 5 or len(args) == 6:
#         pdb_id_one = args[1]
#         pdb_id_two = args[2]
#         anm_mode_number = args[3]
#         chain_id_one = args[4]
#         chain_id_two = args[5] if len(args) == 6 else None
        
#         csv_file_1 = f"{pdb_id_one}_CrossCorr_Mode{anm_mode_number}.csv"
#         csv_file_2 = f"{pdb_id_two}_CrossCorr_Mode{anm_mode_number}.csv"
        
#         csv_path_1 = os.path.join(csv_dir, csv_file_1)
#         csv_path_2 = os.path.join(csv_dir, csv_file_2)
        
#         fasta_path_1 = os.path.join(os.getcwd(), f"{pdb_id_one}.chain{chain_id_one}.nowa.fasta")
#         fasta_path_2 = os.path.join(os.getcwd(), f"{pdb_id_two}.nowa.fasta")
        
#         prot_tile_df_one = prepare_tiles_for_protein(csv_path_1, fasta_path_1)
#         prot_tile_df_two = prepare_tiles_for_protein(csv_path_2, fasta_path_2)
        
#         # Process data between the two proteins
#         result = process_between_proteins_data(prot_tile_df_one, prot_tile_df_two, anm_mode_number, chain_id_one, chain_id_two)
#         # Output the result to CSV
#         result.to_csv(f'filtered_average_cosine_similarities_Mode{anm_mode_number}_between_chain_{chain_id_one}_and_{chain_id_two}.csv', index=False)
        
        
# if __name__ == "__main__":
#     main()