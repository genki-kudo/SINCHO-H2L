import glob
import numpy as np
from matplotlib import pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

def vol_mw_dist(setdir):
    #line up compound's MW and Volume
    comp_list = glob.glob(setdir+'*.pdb')
    #print(comp_list)
    mw_list = []
    vol_list = []
    dist_list = []
    for comp in comp_list:
        #print(comp)
        mol = Chem.MolFromPDBFile(comp, proximityBonding=True)
        mw_list.append(Descriptors.MolWt(mol, onlyHeavy=True))
        vol_list.append(AllChem.ComputeMolVolume(mol, gridSpacing=1.0))
        lst = []
        with open(comp,'r')as inp:
            for line in inp:
                lst.append(np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]))
        dist = []
        for i in range(len(lst)):
            for j in range(len(lst)):
                if i == j:
                    dist.append(0)
                elif i>j:
                    dist.append(float(np.linalg.norm(lst[i] - lst[j])))
                    if float(np.linalg.norm(lst[i] - lst[j]))>=200:
                        print(comp)
                else:pass
                
        dist_list.append(max(dist))

    #linear fitting of Volume-MW
    res1 = np.polyfit(np.array(vol_list), np.array(mw_list), 1)
    res2 = np.polyfit(np.array(dist_list), np.array(mw_list), 1)
    #print(vol_list)
    #print(mw_list)
    #print(dist_list)
    
    
    x = np.linspace(5,10,20)
    y1 = np.poly1d(res2)(x)
    print(y1)
    print(res2[0], res2[1])
    

    plt.scatter(dist_list, mw_list, label='compounds')
    plt.plot(x, y1)
    plt.xlabel('Volume [A^3]')
    plt.ylabel('Molecular Weight')
    plt.legend()
    plt.show()
    
    return res1

