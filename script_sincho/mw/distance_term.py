from script_sincho.basic.basic_func import *

def dist_calc(ligand, poc_and_atom, clusterdir):
    poc = clusterdir+poc_and_atom[0]
    atom = poc_and_atom[1].split('_')[-1]

    for i in open(ligand,'r').readlines():
        if i[12:16].replace(' ','')==atom:
            vec_atom = vec_xyz(i)
    pqr_inf = []
    for i in open(poc,'r').readlines():
        pqr_inf.append(vec_xyz(i))
    
    distance_min = 50000000
    for i in pqr_inf:
        distance_min = dist_cf(vec_atom, i, distance_min)
    
    return distance_min




