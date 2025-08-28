#importing required Libraries
import os
import pandas as pd
import numpy as np
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, NeighborSearch
from Bio import SeqIO

"""Using the downloaded PDB files, we will extract the required information and perform analysis.
Required information and Steps:
1. Parse the PDB files to extract relevant structural information.
2. Remove all waters and heteroatoms.
3. Analyze the distribution of antigen chains and their properties.
4. Select the residues present within 5 Å of the antigen chains and 5 Å of the Heavy and Light chains.
5. Observe the interactions between the antigen chains and the surrounding residues.
6. Observe any network interactions
"""

# pdb_id_path = "D:/VIDO/VIDO PROJECT FILES/PDB_legacy_files/7lz6.pdb"

def parse_pdb(pdb_id_path):
    parser = PDBParser()
    structure = parser.get_structure("PDB_structure", pdb_id_path)
    return structure

def remove_waters_and_heteroatoms(structure):
    for model in structure:
        for chain in model:
            for residue in list(chain):
                if residue.id[0] == "W" or residue.id[0] == "H":
                    chain.detach_child(residue.id)
    return structure

def save_cleaned_pdb(structure, output_path):
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path)

def get_residues_within_distance(structure, chain_id, distance=5.0):
    """Get all residues within a specified distance of a given chain"""
    target_residues = []
    other_residues = []
    
    # Get all atoms from the target chain and other chains
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if chain.id == chain_id:
                        target_residues.append((atom, residue))
                    else:
                        other_residues.append((atom, residue))
    
    # Find residues within distance
    nearby_residues = set()
    
    for target_atom, target_residue in target_residues:
        for other_atom, other_residue in other_residues:
            dist = target_atom - other_atom  # BioPython calculates distance with - operator
            if dist <= distance:
                nearby_residues.add(other_residue)
    
    return list(nearby_residues)

def get_interface_residues_efficient(structure, chain1_id, chain2_id, distance=5.0):
    """More efficient method using NeighborSearch to find interface residues"""
    atoms = []
    chain1_atoms = []
    chain2_atoms = []
    
    # Collect all atoms and separate by chains
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atoms.append(atom)
                    if chain.id == chain1_id:
                        chain1_atoms.append((atom, residue))
                    elif chain.id == chain2_id:
                        chain2_atoms.append((atom, residue))
    
    # Use NeighborSearch for efficient distance calculations
    ns = NeighborSearch(atoms)
    
    interface_residues = set()
    
    # Find chain1 residues near chain2
    for atom, residue in chain1_atoms:
        nearby_atoms = ns.search(atom.coord, distance)
        for nearby_atom in nearby_atoms:
            if nearby_atom.get_parent().get_parent().id == chain2_id:
                interface_residues.add(residue)
                interface_residues.add(nearby_atom.get_parent())
    
    return list(interface_residues)

def get_all_chain_ids(structure):
    """Get all chain IDs from the structure"""
    chain_ids = []
    for model in structure:
        for chain in model:
            chain_ids.append(chain.id)
    return chain_ids

def find_interface_residues_auto(pdb_path, distance=5.0):
    """
    Automatically find interface residues from a cleaned PDB file
    Returns residues within specified distance between different chains
    """
    # Parse the PDB file
    structure = parse_pdb(pdb_path)
    
    # Get all chain IDs
    chain_ids = get_all_chain_ids(structure)
    print(f"Found chains: {chain_ids}")
    
    if len(chain_ids) < 2:
        print("Warning: Need at least 2 chains to find interface residues")
        return []
    
    # Find interface residues between all chain pairs
    all_interface_residues = set()
    
    for i, chain1 in enumerate(chain_ids):
        for chain2 in chain_ids[i+1:]:
            print(f"Finding interface between chains {chain1} and {chain2}")
            interface_residues = get_interface_residues_efficient(structure, chain1, chain2, distance)
            all_interface_residues.update(interface_residues)
            print(f"  Found {len(interface_residues)} interface residues")
    
    print(f"Total unique interface residues: {len(all_interface_residues)}")
    return list(all_interface_residues)

def save_interface_structure(pdb_path, output_path, distance=5.0):
    """
    Extract and save only the interface residues from a PDB file
    """
    # Get interface residues
    interface_residues = find_interface_residues_auto(pdb_path, distance)
    
    if not interface_residues:
        print("No interface residues found")
        return
    
    # Parse original structure
    structure = parse_pdb(pdb_path)
    
    # Create new structure with only interface residues
    interface_residue_ids = set()
    for residue in interface_residues:
        interface_residue_ids.add((residue.get_parent().id, residue.id))
    
    # Remove non-interface residues
    for model in structure:
        for chain in model:
            for residue in list(chain):
                if (chain.id, residue.id) not in interface_residue_ids:
                    chain.detach_child(residue.id)
    
    # Save the interface structure
    save_cleaned_pdb(structure, output_path)
    print(f"Interface structure saved to: {output_path}")

#Extract residues only between 5 Å of the antigen and antibody, only keep the structure that is between them
def extract_residues_within_distance(structure, antigen_chain_id, antibody_chain_ids, distance=5.0):
    """Extract residues within distance between antigen and antibody chains"""
    interface_residues = set()
    
    # Get residues from antibody chains that are within distance of antigen
    antigen_nearby = get_residues_within_distance(structure, antigen_chain_id, distance)
    
    # Get residues from antigen chain that are within distance of antibody chains
    for ab_chain_id in antibody_chain_ids:
        ab_nearby = get_residues_within_distance(structure, ab_chain_id, distance)
        interface_residues.update(ab_nearby)
    
    # Also include antigen residues near antibody
    interface_residues.update(antigen_nearby)
    
    return list(interface_residues)

# test1 = parse_pdb(pdb_id_path)
# test1 = remove_waters_and_heteroatoms(test1)
# save_cleaned_pdb(test1, "D:/VIDO/VIDO PROJECT FILES/pymol_structures/7lz6_cleaned.pdb")


