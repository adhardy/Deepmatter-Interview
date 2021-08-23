from rdkit.Chem.Draw.SimilarityMaps import GetMorganFingerprint
from rdkit.DataStructs import cDataStructs
from rdkit import Chem

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.sparse import csr_matrix
from sklearn.decomposition import TruncatedSVD

from typing import List

from dm_RDFParse import DMMol


def get_molecular_fingerprints(molecules: List[DMMol], nbits: int, sparse_matrix: bool = True): # -> Union[np.array, csr_matrix]:
    """
    Calculates the Morgan fingerprints for a list of molecules.

    == references ===
    convert to numpy - https://iwatobipen.wordpress.com/2019/02/08/convert-fingerprint-to-numpy-array-and-conver-numpy-array-to-fingerprint-rdkit-memorandum/
    rdkit data structure - https://www.rdkit.org/docs/source/rdkit.DataStructs.cDataStructs.html
    sparse matrix - https://machinelearningmastery.com/sparse-matrices-for-machine-learning/
    
    """

    fingerprints = np.zeros((len(molecules),nbits), dtype=np.int8) # intitialize numpy array for fingerprints
    for ii, m in enumerate(molecules):
        try:
            Chem.SanitizeMol(m.mol) # exception in GetMorganFingerprint() without this
            fp = GetMorganFingerprint(m.mol, 2, nBits=1000)

        except:
            print(f"Failed to get fingerprint for molecule {ii}.")

        else:
            arr = np.zeros((0,), dtype=np.int8) # intitialize target numpy array
            cDataStructs.ConvertToNumpyArray(fp, arr)
            fingerprints[ii] = arr
    
    if sparse_matrix:
        fingerprints = csr_matrix(fingerprints) # convert to sparse matrix

    return fingerprints

def get_axes(df: pd.DataFrame, ax: plt.axes, 
    colour_col: str, marker: str = "o", markersize: int = 12, label_idx: bool = True, colors: List[str] = ["tab:blue", "tab:orange"],
    switch_xy: bool = False): # -> plt.figure:

    if not switch_xy: # use if we want to swap the x and y axis
        col_0 = 0
        col_1 = 1
    else:
        col_0 = 1
        col_1 = 0

    groups = df.groupby(colour_col) # group by the molecule type
    for ii, (name, group) in enumerate(groups): # plot each group (reactant/product) - matplotlib will give them different colors
        ax.plot(group[df.columns[col_0]], group[df.columns[col_1]], marker=marker, linestyle='', markersize=markersize, label=name, color = colors[ii])

        ax.legend()
        ax.set_xlabel(df.columns[col_0])
        ax.set_ylabel(df.columns[col_1])

    if label_idx:
    # label each point with it's index in the dataframe/molecules list
        for ii, label in enumerate(df.index):
            ax.annotate(label, 
                    (df.iloc[ii][col_0], df.iloc[ii][col_1]),
                    xytext=(-10, 10),
                    textcoords='offset points')

    return ax

def scree_plot(num_components:int, tSVD: TruncatedSVD, column_names: List[str]):
    fig, ax = plt.subplots(figsize=(num_components/2,5)) # 1 width per PC seems to work about right
    scree = tSVD.explained_variance_ratio_
    ax.bar(x=range(0,num_components), height=scree, tick_label=column_names)
    ax.set_ylabel("Explained Variance (%)")
    ax.set_xlabel("Principal Components")
    return fig