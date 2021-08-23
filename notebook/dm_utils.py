import os
import pathlib

import xml.dom.pulldom as pulldom
import xml.etree.ElementTree as ET

from typing import Any, Callable

import rdkit
from rdkit import Chem

def check_file_extension(file_path: str, extension: str): # -> bool:
    """Checks if a file has a specific extension"""
    if pathlib.Path(file_path).suffix != extension:
        raise TypeError(f"{file_path} is not a .pccr file.")

def none_or_equal(value1:any, value2: Any):# -> bool:
    """Returns True if the values are equal or if either value is None. Otherwise returns False."""
    if value1 is None or value2 is None or value1 == value2:
        return True
    return False

def get_XML_root(pccr_file_path: str, 
            tagname: str, 
            attribute: str = None, 
            attribute_value: Any = None):
    """Generator: yields XML trees matching the given tag and attribute."""

    doc = pulldom.parse(pccr_file_path) # these xml files are very large, use pulldom to extract the parts we need 
    for event, node in doc:
        if event == pulldom.START_ELEMENT and none_or_equal(node.tagName, tagname) and none_or_equal(node.getAttribute(attribute), attribute_value):
            doc.expandNode(node) # expand the node so we can parse it with elementree and xpath
            yield ET.fromstring(node.toxml())  # load the xml into elementree

def process_directory(dir: str, funct: Callable[[str, Any], Any], *args, **kwargs):
    """
    Generator: Runs funct(*args, **kwargs) on all files in the directory.
    
    funct() should have the file_path as the first argument.
    """

    files = os.listdir(dir) # get all files in the directory

    for file in files:
        file_path = os.path.join(dir,file) # get the full path
        yield funct(file_path, *args, **kwargs)

class DMMol():
    """Class to store rdkit mol objects and additional metadata about them."""
    
    def __init__(self, mol:rdkit.Chem.rdchem.Mol, mol_type: str):
        self.mol = mol
        self.mol_type = mol_type

    def smiles(self):
        return Chem.MolToSmiles(self.mol)