import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from pprint import pprint as pp

from scipy.stats import rankdata
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import davies_bouldin_score

import re
import subprocess




def report_tile_indices(file_path, sequence_file_path):
    # Read the CSV file
    path_sequences = pd.read_csv(file_path)

    # Convert Best_Path to boolean if necessary
    path_sequences['Best_Path'] = path_sequences['Best_Path'].astype(bool)

    # Extract the maximum tile length
    max_tile_length = max(path_sequences['Tile_Length'])

    # Define the two points of interest
    half_max_tile_length = 35
    max_tile_length_point = max_tile_length

    # Tolerance for floating-point comparison
    tolerance = 1e-5

    def get_sequence(tile_length, tile_center_index):
        tile_center_index_formatted = f"{tile_center_index:.1f}"
        awk_command = f"gawk -F ',' '{{ if ($5 == \"{tile_center_index_formatted}\" && $9 == \"{tile_length}\") print $1}}' {sequence_file_path}"
        result = subprocess.getoutput(awk_command)
        return result.strip().split("\n")

    def print_sequences(df, tile_length):
        for _, row in df.iterrows():
            sequence = get_sequence(tile_length, row['Tile_Center_Index'])
            print("color gold, ( ps. ", ", ".join(sequence), ")")

    # Filter based on Tile_Length and Best_Path
    half_max_tile_info = path_sequences[
        (abs(path_sequences['Tile_Length'] - half_max_tile_length) < tolerance) & 
        path_sequences['Best_Path']
    ]

    max_tile_info = path_sequences[
        (abs(path_sequences['Tile_Length'] - max_tile_length_point) < tolerance) & 
        path_sequences['Best_Path']
    ]

    print("Half max tile info:")
    print(half_max_tile_info)
    print_sequences(half_max_tile_info, half_max_tile_length)

    print("Max tile info:")
    print(max_tile_info)
    print_sequences(max_tile_info, max_tile_length_point)




def create_sequence_mapping(file_path):
    # Read the CSV file
    df = pd.read_csv(file_path)

    # Create a dictionary mapping (tile_length, tile_center_index) to sequence
    sequence_mapping = {}
    for index, row in df.iterrows():
        tile_length = row['tile_length']
        tile_center_index = row['center_index_1']
        sequence = row['sequence_1']
        sequence_mapping[(tile_length, tile_center_index)] = sequence

    return sequence_mapping



def generate_pymol_script(optimal_sequences, color="gold"):
    pymol_commands = []
    for label, (sequence, path, half_max_length_sequence) in optimal_sequences.items():
        command = f"select ps.{sequence}, color {color}"
        pymol_commands.append(command)
    return "\n".join(pymol_commands)





# def find_optimal_path(island_df, start_col_range=None):
def find_optimal_path(island_df, sequence_mapping, start_col_range=None):
    # Number of rows (Tile_Length) and columns (Tile_Center_Index)
    rows, cols = island_df.shape
    
    # Initialize the DP table
    dp_table = np.full((rows, cols), float('-inf'))
    traceback_table = np.zeros((rows, cols), dtype=int)

    # If start_col_range is provided, initialize only the corresponding columns
    if start_col_range is not None:
        start_col_indices = [island_df.columns.get_loc(col) for col in start_col_range]
        dp_table[0, start_col_indices] = island_df.iloc[0, start_col_indices].values
    else:
        dp_table[0, :] = island_df.iloc[0].values

    # Iterate through the rows (Tile_Length) and columns (Tile_Center_Index)
    for row in range(1, rows):
        for col in range(cols):
            # Values from the previous row and neighboring columns (two on each side)
            previous_values = dp_table[row - 1, max(col - 2, 0):min(col + 3, cols)]
            
            # Maximum value from the previous row's neighboring columns
            max_previous_value = max(previous_values)
            
            # Index of the maximum value (used for traceback)
            max_previous_index = np.argmax(previous_values) + max(col - 2, 0)
            
            # DP recurrence relation
            dp_table[row, col] = island_df.iloc[row, col] + max_previous_value
            
            # Traceback table update
            traceback_table[row, col] = max_previous_index

    # Traceback to find the optimal path
    optimal_path = []
    max_value_index = np.argmax(dp_table[-1])
    optimal_path.append((island_df.index[-1], island_df.columns[max_value_index]))
    for row in range(rows - 1, 0, -1):
        max_value_index = traceback_table[row, max_value_index]
        optimal_path.append((island_df.index[row - 1], island_df.columns[max_value_index]))

    # Reversing the path to start from the smallest Tile_Length
    optimal_path.reverse()

    # Extract the first sequence that satisfies the motif
    optimal_sequence = None
    for tile_length, tile_center_index in reversed(optimal_path):
        sequence = sequence_mapping.get((tile_length, tile_center_index))
        if sequence and validate_motif(sequence):
            optimal_sequence = sequence
            break
        
    # # Determine the index corresponding to half the maximum tile length
    # half_max_length_index = rows // 2
    
    # Determine the index corresponding to half the maximum tile length
    half_max_length_index = len(optimal_path) // 2

    # Extract the sequence at half the maximum length
    half_max_length_tile_length, half_max_length_tile_center_index = optimal_path[half_max_length_index]
    half_max_length_sequence = sequence_mapping.get((half_max_length_tile_length, half_max_length_tile_center_index))

    return optimal_path, optimal_sequence, dp_table[-1, np.argmax(dp_table[-1])], half_max_length_sequence



# Define the motif validation function
def validate_motif(sequence):
    pattern = r'C.*C.*C'
    return re.search(pattern, sequence) is not None


def find_tile_center_index_diff(df):
    # Sort the data by 'Tile_Length' and 'Tile_Center_Index'
    df = df.sort_values(by=['Tile_Length', 'Tile_Center_Index'])

    # Group by 'Tile_Length' and calculate the difference in 'Tile_Center_Index' within each group
    df['Index_Diff'] = df.groupby('Tile_Length')['Tile_Center_Index'].diff().fillna(0.0)

    # Locations where the difference is exactly 2
    locations_with_diff_2 = df[df['Index_Diff'] == 2.0]

    print("\nLocations with Tile Center Index Difference of 2:")
    print("===============================================")

    for tile_length, group in locations_with_diff_2.groupby('Tile_Length'):
        tile_center_indexes = group['Tile_Center_Index'].values
        print(f"Tile Length {tile_length}: {tile_center_indexes}")

    return locations_with_diff_2



def assign_labels(df, merge_ranges_dict, max_tile_length):
    # Unpacking the values from the dictionary
    (pre_tile_length, pre_merge_range1, pre_merge_range2), (post_tile_length, post_merge_range) = next(iter(merge_ranges_dict.items()))

    # Labeling the data based on the pre-merge ranges
    df['Label'] = -1 # Default label

    # Create a mask for the pre-tile length
    pre_tile_mask = df['Tile_Length'] <= pre_tile_length

    # Assign labels for the two pre-merge ranges
    df.loc[pre_tile_mask & df['Tile_Center_Index'].between(*pre_merge_range1), 'Label'] = 0
    df.loc[pre_tile_mask & df['Tile_Center_Index'].between(*pre_merge_range2), 'Label'] = 1

    # Check if there are points still part of label 0 but outside the bounds of the tile right before being merged
    label_0_outside_bounds = pre_tile_mask & (df['Label'] == -1) & (df['Tile_Center_Index'] < pre_merge_range1[0])
    df.loc[label_0_outside_bounds, 'Label'] = 0

    # Check if there are points still part of label 1 but outside the bounds of the tile right before being merged
    label_1_outside_bounds = pre_tile_mask & (df['Label'] == -1) & (df['Tile_Center_Index'] > pre_merge_range1[1]) & (df['Tile_Center_Index'] < pre_merge_range2[0])
    df.loc[label_1_outside_bounds, 'Label'] = 1

    # Labeling the data after the merge
    df.loc[df['Tile_Length'].between(post_tile_length, max_tile_length), 'Label'] = 2

    return df


# Function to create a sequence of points for each path, sorting by Tile_Center_Index and Tile_Length
def create_paths(group):
    group = group.sort_values(by=['Tile_Length', 'Tile_Center_Index'])
    group['Path_ID'] = group['Label'].astype(str) + "_" + group['Starting_Point'].astype(str)
    return group


def identify_and_plot_indices(df):
    boundaries = set([0])  # Start with a boundary of 0

    # Iterate through the tile lengths in the DataFrame
    for tile_length, group in df.groupby('Tile_Length'):
        # Get the tile center indices for the current tile length
        tile_center_indices = group['Tile_Center_Index'].values

        # Iterate through the tile center indices and check for differences
        for i in range(1, len(tile_center_indices) - 1):  # Start from index 1 to check previous difference
            current_difference = tile_center_indices[i + 1] - tile_center_indices[i]
            previous_difference = tile_center_indices[i] - tile_center_indices[i - 1]

            # Check if the current difference is significantly greater than the previous difference
            if current_difference >= 2 * previous_difference and tile_length != 6:
                # Calculate the midpoint between the two indices
                midpoint = (tile_center_indices[i] + tile_center_indices[i + 1]) / 2

                # Plot a vertical line at the midpoint
                plt.axvline(x=midpoint, color='r')

                # Add the midpoint to the boundaries
                boundaries.add(midpoint)

    # Add the maximum Tile_Center_Index value as the last boundary
    boundaries.add(df['Tile_Center_Index'].max() + 1)

    # Convert the set to a sorted list
    boundaries = sorted(list(boundaries))

    # Print the identified boundaries
    print("Boundaries:", boundaries)

    return boundaries



def assign_label_two_nums(value, boundaries):
    for i in range(len(boundaries) - 1):
        if boundaries[i] <= value < boundaries[i + 1]:
            return f'Label {i}'  # Use 'i' instead of 'i + 1'
    return None



def assign_label_two_chars(value, boundaries):
    for i in range(len(boundaries) - 1):
        if boundaries[i] <= value < boundaries[i + 1]:
            return chr(65 + i)  # Returns 'A', 'B', 'C', etc., based on the index
    return None





def main():
    cwd = os.getcwd()
    filename = 'averages_and_counts_ModeAll.csv'
    file_path = os.path.join(cwd, filename)
    
    file_path_seqs = "final_scan_ModeAll_withFalse.csv"
    sequence_mapping = create_sequence_mapping(file_path_seqs)


    df = pd.read_csv(file_path)
    df = df.sort_values(by=['Tile_Length', 'Tile_Center_Index'], ascending=True)

    # * Create a new column to mark non-zero Average_Weighted_Sum
    df['Average_Weighted_Sum_Non_Zero'] = df['Average_Weighted_Sum'] != 0
    non_zero_data_df = df[df['Average_Weighted_Sum_Non_Zero']].copy()

    boundaries = sorted(identify_and_plot_indices(non_zero_data_df))
    non_zero_data_df['Label'] = non_zero_data_df['Tile_Center_Index'].apply(assign_label_two_chars, boundaries=boundaries)

    grid_df = non_zero_data_df.pivot(index='Tile_Length', columns='Tile_Center_Index', values='Average_Weighted_Sum')
    grid_df = grid_df.fillna(float('-inf'))

    print(boundaries)
    print(non_zero_data_df['Label'].unique())


    # Define colormap
    labels = non_zero_data_df['Label'].unique()
    colormap = mpl.colormaps['viridis']
    colormap = colormap(np.linspace(0, 1, len(labels)))
    
    # Create a dictionary to map labels to colors
    label_colors = {label: color for label, color in zip(labels, colormap)}

    # Map the labels to the corresponding colors
    color_values = non_zero_data_df['Label'].map(label_colors)



    # Visualization setup
    plt.figure(figsize=[10, 6])
    plt.scatter(non_zero_data_df['Tile_Center_Index'], non_zero_data_df['Tile_Length'], c=color_values)
    plt.xlabel('Tile Center Index')
    plt.ylabel('Tile Length')
    plt.title('Labels Visualization')

    # Export all labeled points to CSV
    non_zero_data_df.to_csv('all_labeled_points.csv', index=False)

    # Create an empty list to store the path details as individual DataFrames
    path_data_list = []

    # Get unique labels
    unique_labels = non_zero_data_df['Label'].unique()

    # Initialize a dictionary to store the optimal sequences and paths for each label
    optimal_sequences_summary = {}


    # Iterate through each unique label
    for label in unique_labels:
        print(f"\nAnalyzing label: {label}")

        # Filter the data to include only the specified label and the lowest tile length
        label_data_df = non_zero_data_df[non_zero_data_df['Label'] == label]

        # Pivot the data to create a grid where rows represent Tile_Length and columns represent Tile_Center_Index
        grid_df = label_data_df.pivot(index='Tile_Length', columns='Tile_Center_Index', values='Average_Weighted_Sum')
        grid_df = grid_df.fillna(float('-inf'))

        # Determine the starting column range based on tile length 6
        start_col_range = label_data_df[label_data_df['Tile_Length'] == 6]['Tile_Center_Index'].values

        # Determine the optimal sequence and path for the current label
        best_optimal_sequence = None
        best_max_value = float('-inf')
        best_start_col = None
        best_optimal_path = None

        paths_for_label = [] # Initialize a list to store the paths for the current label

        print("Starting Points and Max Values:")
        # Iterate through the starting column range and apply the find_optimal_path function
        for start_col in start_col_range:
            # Apply the find_optimal_path function with the specified starting column
            # optimal_path, max_value = find_optimal_path(grid_df, start_col_range=[start_col])
            # optimal_path, optimal_sequence, max_value = find_optimal_path(grid_df, sequence_mapping, start_col_range=[start_col])
            
            optimal_path, optimal_sequence, max_value, half_max_length_sequence = find_optimal_path(grid_df, sequence_mapping, start_col_range=[start_col])




            if optimal_sequence:
                print("Optimal sequence:", optimal_sequence)
                print(f"Optimal sequence at half the maximum length: {half_max_length_sequence}")
                
            else:
                print("No sequence found that satisfies the motif.")

            print(f"Starting Point: {start_col}, Max Value: {max_value}")  # Print the max value for each starting point

            # Extract Tile_Length and Tile_Center_Index for plotting
            tile_lengths, tile_center_indices = zip(*optimal_path)

            # Plot each optimal path with reduced opacity
            plt.plot(tile_center_indices, tile_lengths, lw=2.5, color="red", alpha=0.3)

            # Create a DataFrame for the current path
            current_path_data = pd.DataFrame({
                'Label': label,
                'Tile_Length': tile_lengths,
                'Tile_Center_Index': tile_center_indices,
                'Starting_Point': start_col,
                'Max_Value': max_value,
                'Best_Path': False, # Initialize with False
                'Path_ID': f"{label}_{start_col}", # Unique path identifier
            })

            # Add an 'order' column that represents the order of points within each path
            current_path_data['order'] = current_path_data.index + 1

            # Add the DataFrame to the list for the current label
            paths_for_label.append((start_col, current_path_data)) # Store the DataFrame along with the starting column


            # # Add the DataFrame to the list for the current label
            # path_data_list.append(current_path_data)

            # Update the best optimal path and maximum value if the current maximum value is greater
            if max_value > best_max_value:
                best_max_value = max_value
                best_optimal_path = optimal_path
                best_start_col = start_col
                best_optimal_sequence = optimal_sequence
                
            # # Print the summary for the current label
            # print(f"Label: {label}")
            # print(f"Best Starting Point: {best_start_col}")
            # print(f"Optimal Sequence: {best_optimal_sequence}")
            # print(f"Maximum Value: {best_max_value}\n")
        
            # # Store the optimal sequence and path for the current label
            # optimal_sequences_summary[label] = (best_optimal_sequence, best_optimal_path)

            # Store the optimal sequence at half the maximum length for the current label
            optimal_sequences_summary[label] = (best_optimal_sequence, best_optimal_path, half_max_length_sequence)


        # Update the 'Best_Path' column for the optimal path
        for start_col, path_df in paths_for_label:
            if start_col == best_start_col:
                path_df['Best_Path'] = True
            path_data_list.append(path_df)  # Append the DataFrame to the path_data_list

            

    # Print a summary of all the labels and their optimal sequences
    print("\nSummary of Optimal Sequences:")
    for label, (sequence, path, half_max_length_sequence) in optimal_sequences_summary.items():
        print(f"Label: {label}")
        print(f"Sequence: {sequence}")
        print(f"Path: {path}")
        print(f"Optimal Sequence at Half the Maximum Length: {half_max_length_sequence}\n")



    # Generate the PyMOL script
    pymol_script = generate_pymol_script(optimal_sequences_summary)
    print(pymol_script)


    # Concatenate all the DataFrames in the list
    path_data = pd.concat(path_data_list, ignore_index=True)

    # Export path_sequences to CSV
    path_data.to_csv('path_sequences.csv', index=False)

    # Invert the y-axis
    plt.gca().invert_yaxis()

    # Show the plot
    plt.show()


    print('\nDone.')


if __name__ == "__main__":
    main()






























