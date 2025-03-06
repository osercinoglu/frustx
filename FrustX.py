import prody as pd 
import sys, os
import subprocess
import argparse
import glob
import tqdm
import re
import panedr
import pandas as pad
import matplotlib.pyplot as plt
import numpy as np
from Bio.SeqUtils import seq1

parser = argparse.ArgumentParser()

pdb_files = glob.glob('proteins/*.pdb')

for i in tqdm.tqdm(range(0,len(pdb_files))):
    pdb_file = pdb_files[i]



    parser.add_argument("--seq")
    parser.add_argument("--chm")
    parser.add_argument("--faspr_path")

    args = parser.parse_args()





    pdb_path = pdb_files[i]
    file_name = pdb_files[i].split("/")
    file_name = file_name[-1].split(".")
    file_place = "proteins/outfolder/"+file_name[0]    
    outFolder = file_place
    seqFile = args.seq
    fasprFile = "FASPR-master/FASPR"

    if args.chm:
        list_a = args.chm
        list_position = list_a.split(",")
        position = 0
    else:
        position = "all"

    firstResidueNo = 0
    myresidue = 0
    mole_org = ""
    # Define the standard amino acids
    amino_acids = 'FLSYCWPHQRIMTNKVADEG'

   
    seqFile = "sequence.txt"

    structure = pd.parsePDB(pdb_path)

    if not os.path.exists(outFolder):
        os.makedirs(outFolder)



    def findFirstNumber(pdb_path):

        # Load the PDB file


        structure = pd.parsePDB(pdb_path)

        #my addition
        residue = structure.select(f'index {0}').getResnums()[0]
        return   residue



    def generate_alternative_sequences(pdb_path, position, firstResidueNo):

        structure = pd.parsePDB(pdb_path)

        #my addition
        native_residue = structure.select(f'resnum {position}').getResnames()[0]
        print("residue number:" , position)

        myresidue = position + 1 - firstResidueNo

        # Assume you want to analyze the first chain (change as needed)
        chain = next(structure.getHierView().iterChains())
        sequence = chain.getSequence()

        # Print the original sequence at the position

        original_residue = sequence[position-firstResidueNo]  # Adjust for zero-indexing
        mole_org = original_residue

        print(f"Original residue at position {position}: {original_residue}")

        # Generate and print alternatives
        alternatives = []
        for aa in amino_acids:
            if aa != original_residue:
                # Create a new sequence with the substitution
                new_sequence = list(sequence)
                new_sequence[position - firstResidueNo] = aa
                new_sequence = ''.join(new_sequence)
                alternatives.append(new_sequence)

        return alternatives, myresidue, mole_org

    if not position == "all":
        for i in list_position:
            position_i = int(i)


            firstResidueNo = findFirstNumber(pdb_path)

            alternative_sequences, myresidue, mole_org = generate_alternative_sequences(pdb_path, position_i, firstResidueNo)

            sequences_list = []

            seq_index = 0

            for seq in alternative_sequences:
                print(seq)
                sequences_list += [seq]
                seq_index += 1




            for seq in sequences_list:
                newFile = mole_org + str(position) + seq[myresidue-1]
                print(os.path.join(os.path.abspath(outFolder),newFile+".pdb"))  
                with open(seqFile, "w") as File:
                    File.writelines(seq)    


                subprocess.run([fasprFile, "-i", pdb_path,"-o", os.path.join(os.path.abspath(outFolder),newFile+".pdb"),"-s", seqFile])
    else:

        structure_CA = structure.select('protein and name CA')
        res_nums = structure_CA.getResnums()
        no_of_residues = len(res_nums)

        for i in range(1,no_of_residues):
            position_i = res_nums[i]
            firstResidueNo = findFirstNumber(pdb_path)
            alternative_sequences, myresidue, mole_org = generate_alternative_sequences(pdb_path, position_i, firstResidueNo)

            sequences_list = []

            seq_index = 0

            for seq in alternative_sequences:
                print(seq)
                sequences_list += [seq]
                seq_index += 1

            for seq in sequences_list:
                newFile = mole_org + str(position_i) + seq[myresidue-1]
                print(os.path.join(os.path.abspath(outFolder),newFile+".pdb"))  
                with open(seqFile, "w") as File:
                    File.writelines(seq)    


                subprocess.run([fasprFile, "-i", pdb_path,"-o", os.path.join(os.path.abspath(outFolder),newFile+".pdb"),"-s", seqFile])

    WT_list = []
    Wild_Type = "WT_file.txt"
    # Load the PDB file

    pdb_WT = pdb_file
    seq_2 = pd.parsePDB(pdb_WT)
    
    # Select only standard amino acid residues (protein chains)
    protein = seq_2.select('protein and name CA')  # Selects only alpha carbon atoms (removes duplicates)

    # Extract residue names and convert them to one-letter codes
    if protein:
        residues = protein.getResnames()
        sequence = ''.join(seq1(res) for res in residues)
        with open(Wild_Type, "w") as File:
            File.writelines(sequence)       
        print(f"Corrected Protein Sequence: {sequence}")
    else:
        print("No protein residues found in the PDB file.")
        
    subprocess.run([fasprFile,"-i", pdb_WT, "-o", os.path.join(os.path.abspath(outFolder),"WT.pdb"), "-s", Wild_Type])


                
# List all pdb files in all subfolders in single_aa_variants parent folder.
pdb_files = glob.glob('proteins/*/*/*.pdb')

for i in tqdm.tqdm(range(0,len(pdb_files))):
    pdb_file = pdb_files[i]
    #if "WT" not in pdb_file:
        #continue
    
    pdb_id = re.search(r'\/(\w+).pdb', pdb_file).group(1)
    folder = file_name[0] 

    # Create outfolder in grinn_workflow_output folder.
    outfolder = f'grinn_workflow_output/{folder}/{pdb_id}'

    # Run grinn_workflow.py for each pdb file.
    p = subprocess.run(['python', 'repos/grinn/grinn_workflow.py', pdb_file, 'repos/grinn/mdp_files', outfolder, '--nointeraction', '--gmxrc_path','/opt/gromacs2023_4/bin/GMXRC', '--nt', '24'])
    
    
# Parse all *.edr files in all subfolders of grinn_workflow_output
edr_files = glob.glob('grinn_workflow_output/**/**/*.edr', recursive=True)

# Parse each file to get a pandas dataframe
dfs = []
for edr_file in tqdm.tqdm(edr_files):
    try:
        df = panedr.edr_to_df(edr_file)

        # Get the folder name of the file, which includes mutation in the format: "oldres_position_newres" in aa one-letter codes

        folder_list = edr_file.split('/')
        
        folder = folder_list[1]+'/'+folder_list[2]
        pdb = folder_list[1]

        if "WT" in folder:
            df['mutation'] = 'WT'
            df['oldres'] = 'WT'
            df['newres'] = 'WT'
            df['position'] = 'WT'
        else:
            mutation = re.search(r'[A-Z][0-9]+[A-Z]', folder).group(0)
            df['mutation'] = mutation
            # Further split mutation into parts to get oldres, position, and newres in individual variables
            oldres = mutation[0]
            position = int(mutation[1:-1])
            newres = mutation[-1]
            df['oldres'] = oldres
            df['position'] = position
            df['newres'] = newres
        
        df['pdb'] = pdb
        df['folder'] = folder

        # Take only the last row
        df.drop(index=df.index[:-1],axis=0, inplace=True)
        dfs.append(df)
    except:
        print(f"Failed to parse {edr_file}")

# Concatenate all dataframes into a single dataframe
if dfs:
    df_merged = pad.concat(dfs, axis=0)
else:
    print("No valid .edr files were processed.")
    exit()

# Save the dataframe to a csv file
df_merged.to_csv('grinn_output_table.csv', index=False)

df_stripped = df_merged[['Potential', 'pdb', 'oldres', 'position', 'newres', 'folder', 'mutation']]
df_stripped.to_csv('grinn_stripped.csv', index=False)

# Load CSV
df = pad.read_csv("grinn_stripped.csv", sep=",")

# Display all columns and rows
pad.set_option("display.max_columns", None)
pad.set_option("display.max_rows", None)

# Get unique positions
uni_value = df["position"].unique().tolist()

# Function to calculate the average potential at a given position
def findAverage(i):
    temp_df = df[df["position"] == i]
    return temp_df["Potential"].mean()

# Function to get WT Potential for a given PDB value
def findWT(pdb_value):
    temp_pdb = df[df["pdb"] == pdb_value]
    wtLine = temp_pdb[temp_pdb["position"] == "WT"]
    
    if not wtLine.empty:
        return wtLine["Potential"].values[0]
    return np.nan  # Return NaN if 'WT' not found

# Initialize frustration_index column with NaN
df["frustration_index"] = np.nan

# Iterate over each unique position
for i in uni_value:
    N = df["position"].value_counts().get(i, 0)
    residueDF = df[df["position"] == i]
    
    mySigma = 0  # Reset mySigma for each residue position
    
    # First loop: Calculate variance component (σ²)
    for index, row in residueDF.iterrows():
        WT = findWT(row["pdb"])
        if np.isnan(WT):  # Skip if WT value is missing
            continue
        mySigma += (row["Potential"] - findAverage(i)) ** 2

    # Second loop: Calculate frustration index
    for index, row in residueDF.iterrows():
        WT = findWT(row["pdb"])
        if np.isnan(WT) or mySigma == 0 or N == 0:  # Avoid division by zero
            continue
        result1 = (findAverage(i) - WT) / (np.sqrt(mySigma / N))
        df.loc[row.name, "frustration_index"] = result1  # Correctly update main DataFrame
        
df.to_csv("indexes.csv", sep=",", index=False)  # Saves without the index column