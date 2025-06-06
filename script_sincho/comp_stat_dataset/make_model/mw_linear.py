
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors, AllChem
import plotly.graph_objects as go
import numpy as np
from sklearn import datasets
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
import seaborn as sns
import pandas as pd
from rdkit.Chem.Crippen import MolLogP

def visual_2d(xlist, ylist, xlabel, ylabel, output, diagonal=False):
    fig = plt.figure(figsize=(20,16))
    plt.scatter(xlist, ylist, s=10)
    correlation_coefficient = np.corrcoef(xlist, ylist)[0, 1]
    print(correlation_coefficient)
    plt.xlabel(xlabel,fontsize=35)
    plt.ylabel(ylabel,fontsize=35)
    plt.tick_params(labelsize=30)
    if diagonal:
        plt.plot([min(xlist),max(xlist)],[min(xlist),max(xlist)], color="gray", linestyle="dashed")
    plt.savefig(output)

def est_mw_ind(comp,spath, sbl_list):
    #spath = Chem.GetShortestPath(mol, 5, 10)
    
    atom = comp.GetAtoms()
    print(atom)
    est_mw_ind = 0
    for a in atom:
        if a.GetIdx() in spath[1:-1]:
            #print(a.GetMass())
            sbl_list.append(a.GetSymbol())
            est_mw_ind+=a.GetMass()
    return(est_mw_ind)

#dfをcsvから読み込み
df = pd.read_csv("descriptors4.csv")
print(df.shape)


for index, line in enumerate(df.iterrows()):
    id = line[1][4]
    #print(id)
    mol = Chem.SDMolSupplier("../20241202_PDBbind_v2020_refined/refined-set/"+id+"/"+id+"_ligand.sdf",removeHs=False, sanitize=False)[0]
    """
    if Chem.rdMolDescriptors.CalcNumLipinskiHBA(mol)>=10 or Chem.rdMolDescriptors.CalcNumLipinskiHBD(mol)>=5:
        #or Chem.Descriptors.MolWt(mol)>=500:
        #or float(line[1][22])/float(line[1][23])>=0.5:
        df.drop( df [ df["PDB"].str.contains(id) ] .index, inplace=True )
    """

print(df.shape)

"""
pdbids = []
distances = []
passmws = []
volume = []
allmw = []


for index, line in enumerate(df.iterrows()):
    id = line[1][3]
    print(id)
    sdf = "../20241202_PDBbind_v2020_refined/refined-set/"+id+"/"+id+"_ligand.sdf"
    mol = Chem.SDMolSupplier(sdf,removeHs=False, sanitize=False)[0]

    conformer = mol.GetConformer()
    num_atoms = mol.GetNumAtoms()
    dist = 0
    atom1, atom2 = -1, -1

    for i in range(num_atoms):
        for j in range(i+1, num_atoms):
            pos1 = np.array(conformer.GetAtomPosition(i))
            pos2 = np.array(conformer.GetAtomPosition(j))
            dist_tmp = np.linalg.norm(pos1 - pos2)
            if dist_tmp >= dist:
                dist = dist_tmp
                atom1, atom2 = i,j



    lists = [dist, atom1, atom2]
    print(lists)

    spath = Chem.GetShortestPath(mol, atom1, atom2)
    pass_mw = est_mw_ind(mol,spath, [])
    
    pdbids.append(id)
    distances.append(dist)
    passmws.append(pass_mw)
    volume.append(AllChem.ComputeMolVolume(mol, gridSpacing=1.0))
    allmw.append(Descriptors.MolWt(mol, onlyHeavy=True))
"""

pdbids = df["PDB"].tolist()
distances = df["dist"].tolist()
passmws = df["pass_mw"].tolist()
volume = df["volume"].tolist()
allmw = df["mw"].tolist()

visual_2d(distances, passmws, "Max-Distance in Cpd. [Å]", "MW of Shortest Path in Cpd.", "dist-mw.png", diagonal=False)
res1 = np.polyfit(np.array(distances), np.array(passmws),1)
corr_dist_passmw = np.corrcoef(np.array(distances), np.array(passmws))[0, 1]
print("dist-mw", res1)
print("Correlation coefficient (dist-passmw):", corr_dist_passmw)

visual_2d(volume, allmw, "Volume of Cpd. [Å$^{3}$]", "MW of Cpd.", "vol-mw.png", diagonal=False)
res1 = np.polyfit(np.array(volume), np.array(allmw),1)
corr_vol_mw = np.corrcoef(np.array(volume), np.array(allmw))[0, 1]
print("Correlation coefficient (vol-mw):", corr_vol_mw)
print("vol-mw", res1)

