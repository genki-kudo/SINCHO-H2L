import os
import itertools
from script_sincho.basic.basic_func import *
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

"""
"""
#test-1 4YK0_traj10
cdir = "/home/user01/Research/00_AIH2L/extend-point-selection/\
230710_ensembles_treatment/\
4YK0/md/separate_file/traj_10/"
poc = cdir+"p2c_output/cluster/cluster3.pqr"
lig = cdir+"lig_10.pdb"
lead = "./leadst.pdb"
anchoratom = "C12"
#poc_lat
lat_gen(poc, "ATOM  ", "./lat_.pdb")
latveclist = []
listing_pdbcoordinate("./lat_.pdb", latveclist)
anchor = coordinate_of_atom(lig, "C7")

"""
#test-2 2UW3_crystal
cdir = "/home/user01/Research/00_AIH2L/extend-point-selection/\
230710_ensembles_treatment/\
2UW3/"
poc = cdir+"p2c_output/cluster/cluster6.pqr"
lig = cdir+"HITh.pdb"
lead = "./2uw7.pdb"
anchoratom = "C14"
#poc_lat
lat_gen(poc, "ATOM  ", "./lat_2.pdb")
latveclist = []
listing_pdbcoordinate("./lat_2.pdb", latveclist)
anchor = coordinate_of_atom(lig, "C11")

#test-3 6AXQ_crystal
cdir = "/home/user01/Research/00_AIH2L/extend-point-selection/\
230710_ensembles_treatment/\
6AXQ/"
poc = cdir+"p2c_output/cluster/cluster6.pqr"
lig = cdir+"HITh.pdb"
lead = "./6AY3_l.pdb"
anchoratom = "C12"
#poc_lat
lat_gen(poc, "ATOM  ", "./lat_3.pdb")
latveclist = []
listing_pdbcoordinate("./lat_3.pdb", latveclist)
anchor = coordinate_of_atom(lig, "C7")

#test-4 5FNQ_1_crystal
cdir = "/home/user01/Research/00_AIH2L/extend-point-selection/\
230710_ensembles_treatment/\
5FNQ/"
#poc = cdir+"p2c_output/cluster/cluster11.pqr"
poc = cdir+"p2c_output/cluster/cluster9.pqr"
lig = cdir+"HITh.pdb"
#lead = "./5FNU_1.pdb"
lead = "./5FNU_2.pdb"
#anchoratom = "C13"
anchoratom = "C23"
#poc_lat
#lat_gen(poc, "ATOM  ", "./lat_4.pdb")
lat_gen(poc, "ATOM  ", "./lat_5.pdb")
latveclist = []
#listing_pdbcoordinate("./lat_4.pdb", latveclist)
listing_pdbcoordinate("./lat_5.pdb", latveclist)
#anchor = coordinate_of_atom(lig, "C5")
anchor = coordinate_of_atom(lig, "C11")

"""


#poc's axes
main_axes_val = 0
for l in latveclist:
    distance = float(np.linalg.norm(l - anchor))
    if distance >= main_axes_val:
        main_axes_val = distance
        rep_atom = l
main_axes_vec = rep_atom-anchor

sub_axes_val =0
for pair in itertools.combinations(latveclist,2):
    p = np.array(list(pair)[0]) - np.array(list(pair)[1])
    distance = float(np.linalg.norm(p))
    x = np.inner(main_axes_vec, p)
    a = np.linalg.norm(main_axes_vec)
    b = np.linalg.norm(p)
    theta = np.arccos(x/(a*b))
    dis = distance * np.sin(theta)
    if dis>=sub_axes_val:
        sub_axes_val = dis
        t = theta
        sub_axes_vec = list(pair)

#print(main_axes_vec, main_axes_val)
#print(sub_axes_vec, sub_axes_val)
print("anchor_pocdist->", main_axes_val)
print("ortho_pocdist->", sub_axes_val)
print("I3->", abs(main_axes_val-sub_axes_val)/min([main_axes_val,sub_axes_val]))
print("e->", np.sqrt(1-((0.5*min([main_axes_val,sub_axes_val]))**2/(0.5*max([main_axes_val,sub_axes_val]))**2)))

print("")
#lig's axes
idx = 0
for f in open(lead):
    if f[0:6]=="ATOM  "or f[0:6]=="HETATM":
        if f[12:16].replace(' ','') == anchoratom:
            anchor_idx=idx
        idx+=1

mol = Chem.MolFromPDBFile(lead, removeHs=True)
atoms = mol.GetNumAtoms()

mains = max([len(Chem.GetShortestPath(mol, anchor_idx, i)) for i in range(atoms) if i != anchor_idx])

print("anchor_maxpath->", mains)

l=[]
for i in range(atoms):
    l.append(max([len(Chem.GetShortestPath(mol, i, j)) for j in range(atoms) if j!=i]))
subs = min(l)
print("min_maxpath->",subs)
print("I3->", (mains-subs)/subs)
print("e->", np.sqrt(1-((0.5*subs)**2/(0.5*mains)**2)))
print("")



"""
for i in range(atoms):
    for j in 

idx = mol.GetAtoms()
print(idx)
for i in idx:
    print(i)

"""








