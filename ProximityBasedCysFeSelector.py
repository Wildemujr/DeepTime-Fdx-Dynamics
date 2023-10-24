import os
import sys
import pandas as pd
import numpy as np
from collections import OrderedDict
from Bio.PDB import PDBParser



def calculate_cysteine_positions(row):
    """
    Calculate the positions of cysteines in a given sequence.

    Args:
        row (pd.Series): A row from the DataFrame containing sequence data.

    Returns:
        list: Positions of cysteines in the sequence.
    """
    sequence = row['sequence_1']
    start_index = row['start_index_1']
    positions = [i + start_index for i, x in enumerate(sequence) if x == 'C']
    return positions


def calc_distance(atom_1, atom_2):
    """
    Calculate the distance between two atoms.

    Args:
        atom_1, atom_2 (np.array): The coordinates of the two atoms.

    Returns:
        float: The distance between the two atoms.
    """
    diff_vector = atom_1 - atom_2
    return np.round(np.sqrt(np.sum(diff_vector**2)), decimals=2)


def write_pml_file(pdb, chain_id, fname, mode, directory):
    """
    Write a .pml file based on the DataFrame and mode.

    Args:
        pdb (str): The pdb name.
        chain_id (str): The chain id.
        fname (str): The filename to write to.
        mode (str): The mode.
        directory (str): The current working directory.
    """
    df = pd.read_csv(f"final_scan_Mode{mode}.csv")
    df = df.drop_duplicates(subset='sequence_1', keep='first')

    with open(fname, "w+") as f:
        f.write(f"load {directory}/{pdb}.chain{chain_id}.nowa.pdb\n")
        f.write(f"color gray90, {pdb}.chain{chain_id}.nowa\n")
        
        for sequence in df['sequence_1']:
            f.write(f"select m. {pdb}.chain{chain_id}.nowa and ps. {sequence}\n")
            f.write("color gold, sele\n")




def main():
    # Parse command line arguments for mode number, structure ID, and chain ID
    mode_num = sys.argv[1]
    struct_id = sys.argv[2]
    chain_id = sys.argv[3]
    
    # Get the current working directory
    directory = os.getcwd()

    # Load the cosine similarity data into a DataFrame
    df = pd.read_csv(f"cosine_similarities_Mode{mode_num}_updated_ranked_LocalAln_weightedSum_whole.csv")

    # Check sequences for the presence of cysteines ('C')
    df['Contains_C'] = df['sequence_1'].str.contains("C")
    df['Contains_C'] = df['Contains_C'].replace({True: 'TRUE', False: 'FALSE'})
    
    # Calculate the positions of cysteines in each sequence
    df["Cys_Positions"] = df.apply(calculate_cysteine_positions, axis=1)

    # Initialize PDB parser and load structures
    parser = PDBParser()
    structure = parser.get_structure(struct_id, f"{struct_id}.chain{chain_id}.nowa.pdb")
    hetatm_structure = parser.get_structure(f"{struct_id}_hetatm", f"{struct_id}.chain{chain_id}.nowa.Cofactors.pdb")

    # Extract iron atoms from the bacterial iron-sulfur cluster residue (SF4)
    iron_atoms = [atom for atom in hetatm_structure.get_atoms() if atom.get_parent().get_resname() == 'SF4' and atom.get_name() in [f'FE{num}' for num in range(1, 5)]]

    # Identify permissible cysteines based on their distance to iron atoms
    permissible_cysteines = [
        res.id[1] for res in structure.get_residues() if res.resname == "CYS" and any(
            2.1 <= calc_distance(atom.get_coord(), iron_atom.get_coord()) <= 2.5 
            for atom in res if atom.get_name() == "SG" 
            for iron_atom in iron_atoms
        )
    ]

    # Remove duplicates from the permissible cysteines list
    permissible_cysteines = sorted(list(OrderedDict.fromkeys(permissible_cysteines)))
    
    # Mark sequences with permissible cysteine positions for retention
    df['Positions_to_Keep'] = df['Cys_Positions'].apply(lambda x: any(elem in permissible_cysteines for elem in x))
    df['Positions_to_Keep'] = df['Positions_to_Keep'].replace({True: 'TRUE', False: 'FALSE'})
    
    # Set the weighted sum to zero for sequences without permissible cysteines
    df.loc[df['Positions_to_Keep'] == "FALSE", 'weighted_sum'] = 0
    df.to_csv(f"final_scan_Mode{mode_num}_withFalse.csv", index=False)

    # Filter sequences to only those with permissible cysteines
    df = df[df['Positions_to_Keep'] == "TRUE"]
    if df.empty:
        print(f"final_scan_Mode{mode_num}.csv is empty.")
    else:
        df.to_csv(f"final_scan_Mode{mode_num}.csv", index=False)
        
        # Generate the .pml file for visualization
        write_pml_file(struct_id, chain_id, f"filtered_scan_forCys_Mode{mode_num}.pml", mode_num, directory)

    
if __name__ == '__main__':
    main()
