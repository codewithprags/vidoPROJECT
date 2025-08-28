"""The goal of the PDB download script is to automate the process of downloading PDB legacy files from the RCSB PDB website based on a list of PDB IDs.
Main GOALS:
1. Download PDB legacy files for a list of PDB IDs.
2. Save the downloaded PDB legacy files to a specified directory.
3. Handle errors and retries for failed downloads.
4. Provide progress feedback to the user.
5. Ensure the files are stored as .pdb files that can be easily edited to remove excess information and chains"""


#download required libraries###

from Bio import PDB
from Bio.PDB import PDBParser, PDBIO
from Bio import SeqIO
import numpy as np
import pandas as pd
import requests
import matplotlib.pyplot as plt
import seaborn as sns

def subset_pdb_ids(file_path):
    """
    Select PDB IDs from the main csv file , select first 2000 PDB IDS as a csv file
    """
    selected_ids = []
    df = pd.read_csv(file_path)
    for index, row in df.iterrows():
        if index < 2000:
            selected_ids.append(row['pdb'])

    subset_df = df[df['pdb'].isin(selected_ids)]
    subset_df.to_csv("subset_pdb_ids.csv", index=False)



    #return the csv file containing all values for the selected pdb IDS
    return "subset_pdb_ids.csv"


#download the PDB legacy files for the selected PDB ids
def download_pdb_files(pdb_ids, save_directory):
    for pdb_id in pdb_ids:
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(url)
        if response.status_code == 200:
            with open(f"{save_directory}/{pdb_id}.pdb", "wb") as f:
                f.write(response.content)
            print(f"Downloaded {pdb_id}")
        else:
            print(f"Failed to download {pdb_id}")


#Select the PDB IDs from the subset csv file , filter for NAN, MISSING and  multiple antigen chain.
#check each pdb unique ID, if there are multiple antigen chains, keep only the first one that does not have multiple antigen chains separated by "|"

def filter_pdb_ids(file_path):
    df = pd.read_csv(file_path)
    
    print(f"Original rows: {len(df)}")
    
    # Remove rows with missing values in required columns
    filterdf = df.dropna(subset=['pdb', 'antigen_chain', 'Hchain', 'Lchain'])
    print(f"After removing NaN values: {len(filterdf)}")
    
    # Remove rows where antigen_chain contains "|"
    filterdf = filterdf[~filterdf['antigen_chain'].str.contains('|', regex=False, na=False)]
    print(f"After removing rows with '|' in antigen_chain: {len(filterdf)}")
    
    # Keep only the first occurrence of each unique pdb id
    filterdf = filterdf.drop_duplicates(subset=['pdb'], keep='first')
    print(f"After keeping first occurrence of each PDB: {len(filterdf)}")
    
    # Save the filtered data
    filterdf.to_csv("filtered_pdb_ids.csv", index=False)
    
    print(f"Saved {len(filterdf)} unique PDB entries to filtered_pdb_ids.csv")


    return "filtered_pdb_ids.csv"