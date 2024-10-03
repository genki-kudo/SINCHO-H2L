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

def polar_ratio(pdbfile):
    fullpolar = ["ASP","GLU","ARG","LYS","HIS","HIE","HID","HIP"]
    halfpolar = ["THR","SER","ASN","GLN","TYR","CYS","CYX","XXX"]

    mc = ["CA ","C  ","N  ","O  "]
    alist = []
    count, p_count = 0,0
    for atom in open(pdbfile):
        a = atom[13:26]
        if a[0:3] in mc:
            alist.append("GLY"+a[7:])
        else:
            alist.append(a[4:])
    for r in alist:
        count+=1
        if r[0:3] in fullpolar:
            p_count+=1
        elif r[0:3] in halfpolar:
            p_count+=0.5
    estimate_pocenv = str(p_count/count)

    return estimate_pocenv




