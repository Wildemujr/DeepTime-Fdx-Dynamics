import pandas as pd
import os
import sys


def load_data(anm_mode_number, directory):
    """
    Load data from the CSV file.

    Args:
        anm_mode_number (str): The mode number for the ANM.
        directory (str): The current working directory.

    Returns:
        pd.DataFrame: Loaded data.
    """
    csv_file_path = os.path.join(directory, f"final_scan_Mode{anm_mode_number}_withFalse.csv")
    return pd.read_csv(csv_file_path)


def calculate_averages_and_counts(df):
    """
    Calculate average values and counts for the specified columns.

    Args:
        df (pd.DataFrame): Input data.

    Returns:
        pd.DataFrame: Data with averages and counts.
    """
    # Group by 'center_index_1' and 'tile_length'
    grouped = df.groupby(['center_index_1', 'tile_length'])
    
    # Calculate averages and counts
    averages = grouped.agg({
        'cosine_similarity': 'mean',
        'sequence_similarity': 'mean',
        'frobenius_similarity': 'mean',
        'frobenius_similarity_normed': 'mean',
        'weighted_sum': 'mean',
        'center_index_1': 'size'
    })
    
    averages.rename(columns={
        'cosine_similarity': 'Average_Cosine_Similarity',
        'sequence_similarity': 'Average_Sequence_Similarity',
        'frobenius_similarity': 'Average_Frobenius_Similarity',
        'frobenius_similarity_normed': 'Average_Frobenius_Similarity_Normed',
        'weighted_sum': 'Average_Weighted_Sum',
        'center_index_1': 'Count'
    }, inplace=True)
    
    return averages.reset_index()


def main():
    # Ensure proper command-line arguments are provided
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <anm_mode_number>")
        sys.exit(1)

    anm_mode_number = sys.argv[1]
    directory = os.getcwd()

    df = load_data(anm_mode_number, directory)
    result = calculate_averages_and_counts(df)

    # Rename columns for clarity
    result.rename(columns={
        'center_index_1': 'Tile_Center_Index',
        'tile_length': 'Tile_Length'
    }, inplace=True)

    # Save the result to a CSV file
    result.to_csv(f'averages_and_counts_Mode{anm_mode_number}.csv', index=False)


if __name__ == '__main__':
    main()
