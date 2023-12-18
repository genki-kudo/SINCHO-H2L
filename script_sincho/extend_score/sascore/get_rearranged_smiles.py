import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
import os
import glob
import sys
import rdkit
from rdkit import Chem
#from IPython.core.debugger import Pdb
from openbabel import pybel
#from openbabel import pybel
import openbabel

print(openbabel.__version__)


def rearrange_smiles(smi,atom_idx):
    pbmol = pybel.readstring('smi', smi)
    conv = pybel.ob.OBConversion()
    conv.SetOutFormat("smi")
    conv.SetOptions('l"%d"'%(atom_idx), conv.OUTOPTIONS)     # 1始まりなので+1
    rearranged_smiles = conv.WriteString(pbmol.OBMol).split()[0]
    # print(f"({atom_idx) " + rearranged_smiles)   # 出力文字列の最後に"\t\n"が付いていたのでsplitで切り離し
    return rearranged_smiles

def read_mol(pdb_path):
    mol1=Chem.MolFromPDBFile(pdb_path,sanitize=False)
    # smi = Chem.MolToSmiles(mol,isomericSmiles=True)
    smi = Chem.MolToSmiles(Chem.rdmolops.RemoveHs(mol1),isomericSmiles=True)
    mol2 = Chem.MolFromSmiles(smi)
    smi = Chem.MolToSmiles(Chem.rdmolops.RemoveHs(mol2),isomericSmiles=True)
    for mol in [mol1,mol2]:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx()+1)
    mol1=Chem.rdmolops.RemoveHs(mol1)

    map_num = []
    for atom in mol1.GetAtoms():
        map_num.append(atom.GetAtomMapNum())

    mat = list(mol2.GetSubstructMatch(mol1))
    mat = [m+1 for m in mat]
    match = pd.DataFrame(mat).reset_index(drop=False)
    match['index'] = match['index'] +1
    match['H_num_index'] = map_num
    match.columns = ['no_H_num_index','obabel_num','H_num_index']
    match = match[['obabel_num','no_H_num_index','H_num_index']]

    return match,smi

def get_obabel_num(match,extend_idx):
    extend_idx = int(extend_idx)
    return int(match[match['H_num_index']==extend_idx]['obabel_num']) 

if __name__ == "__main__":
    args = sys.argv
    extend_idx = args[1]
    pdb_path = args[2]
    mol_dir = os.path.dirname(pdb_path)

    match,smi = read_mol(pdb_path)

    obabel_num = get_obabel_num(match,extend_idx)

    smiles = rearrange_smiles(smi,obabel_num)

    print(smi)

    recompound_file = os.path.join(mol_dir,'recompound.txt')
    with open(recompound_file,mode='w') as f:
        f.write(smiles)
