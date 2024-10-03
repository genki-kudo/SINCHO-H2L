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

def around_as(pqrfile, protein, distance, outputdir):

    aminoacidlist =["GLY","ALA","SER","THR","VAL","ILE","LEU","MET","CYS","CYX","HID",'HIE',"HIP"
                   "PHE","HIS","TYR","TRP","PRO","ARG","LYS","ASN","ASP","GLN","GLU"]

    as_dist_file = outputdir+"/pocket_environment/pocenv"+pqrfile.split("/")[-1].split(".")[0][7:]+".pdb"
    t_file(as_dist_file)
    as_xyz = [vec_xyz(i) for i in open(pqrfile) if len(vec_xyz(j))==3]
    for j in open(protein):
        if (j[0:6]=="ATOM  " or j[0:6]=="HETATM") and j[17:20] in aminoacidlist:
            mindist = 1000000
            for a in as_xyz:
                dist = np.linalg.norm(a-np.array([float(j[30:38]),float(j[38:46]),float(j[46:54])]))
                if dist<= mindist:
                    mindist = dist
                if mindist<=distance:
                    with open(as_dist_file,'a')as d:
                        print(j, end="", file=d)
                    break
    return as_dist_file





