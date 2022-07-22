import glob
import numpy as np
from matplotlib import pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

def vol_mw(setdir):
    #line up compound's MW and Volume
    comp_list = glob.glob(setdir+'*.pdb')
    mw_list = []
    vol_list = []
    for comp in comp_list:
        mol = Chem.MolFromPDBFile(comp, proximityBonding=True)
        mw_list.append(Descriptors.MolWt(mol))
        vol_list.append(AllChem.ComputeMolVolume(mol, gridSpacing=1.0))
    #linear fitting of Volume-MW
    res1 = np.polyfit(vol_list, mw_list, 1)
    #print(vol_list)
    #print(mw_list)
    """
    x = np.linspace(100,500,800)
    y1 = np.poly1d(res1)(x)
    print(y1)

    plt.scatter(vol_list, mw_list, label='compounds')
    plt.plot(x, y1)
    plt.xlabel('Volume [A^3]')
    plt.ylabel('Molecular Weight')
    plt.legend()
    plt.show()
    """
    return res1

