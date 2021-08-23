from typing import List


from rdkit.Chem import rdChemReactions
from rdkit import Chem

from dm_utils import DMMol

def extract_reactions(file_path:str): # -> List[str]:
    """Extracts the reactions from a .rdf file and returns them as a list of strings."""

    section_seperator = "$RXN"
    rxn_blocks = [] # this will store each extracted reaction block
    rxn_block = "" # temporarily storos individual reaction blocks while they are parsed
    in_header = True # indicates if we are in the header of the file - we ignore the header
    eof = False # indicates if we are at the end of the file

    # === Parse File ===
    with open(file_path, "r") as f:
        while eof == False:
            line = f.readline()

            if line == "":
                eof = True

            if line == f"{section_seperator}\n":
                if in_header == True: 
                    in_header = False # if we are here, we have found the first reaction block
                else:
                    rxn_blocks.append(rxn_block) # append the previous block to the list of reactions
                    rxn_block = "" # start a new block
            
            if in_header == False:
                rxn_block += line

    return rxn_blocks

def get_molecules_from_rxn(rxn:str):
    """Extracts the reactants and products from a $RXN string."""

    # rdChemReactions - https://www.rdkit.org/docs/source/rdkit.Chem.rdChemReactions.html
    products = rxn.GetProducts() # get all products
    reactants = rxn.GetReactants() # get all reactants
    molecules = []

    for reagents, reagent_type in zip([reactants, products], ["reactant", "product"]):
        # put all reactants and products in one list, using the DMMol class to store the mol object and the type of reagent (reactant or product)
        molecules += [DMMol(reagent, reagent_type) for reagent in reagents]

    return molecules

def get_molecules_from_rxn_list(rxn_blocks:list):
    """Extracts the reactants and products from a list of $RXN strings."""

    max_reactions = 5
    molecules = [] # list of DMMol objects

    for ii in range(max_reactions):
        
        rxn = rdChemReactions.ReactionFromRxnBlock(rxn_blocks[ii]) # parse RXN block string using rdkit
        rxn_molecules = get_molecules_from_rxn(rxn) # get all reactants and products from that reaction
        molecules += rxn_molecules
    
    return molecules

def print_molecules_from_molecule_list(molecules: List[str]):
    
    for molecule in molecules:
        print(f"{molecule.mol_type}: {molecule.smiles()}") # print SMILES string of each molecule
