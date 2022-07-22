from script.basic_func import t_file
from script.basic_func import vec_xyz
import numpy as np


def r_candidate(ligand, candidates):
    t_file(candidates)
    atomnumber_list = []
    with open(ligand,'r')as lig1:
        for line in lig1:
            if (line[0:6]=="HETATM" or line[0:6]=="ATOM  ") and line[76:78]!=' H':
                atomnumber_list.append(int(line[7:11]))
    #print(atomnumber_list)
    candidate =[]
    atom_hydro_vec ={}
    with open(ligand,'r')as lig1:
        for line in lig1:
            if line[0:6]=="CONECT" and int(line[7:11]) not in atomnumber_list:
                line_split = line.split()
                candidate.append(line_split[2])
                hydro_vec = []
                ha_vec = []
                with open(ligand, 'r')as lig2:
                    for line2 in lig2:
                        if line2[0:6]=="HETATM" or line2[0:6]=="ATOM  ":
                            if str(int(line2[7:11]))==str(int(line_split[1])):
                                hydro_vec = vec_xyz(line2)
                            elif str(int(line2[7:11]))==str(int(line_split[2])):
                                ha_label = line2[13:16]
                                ha_vec = vec_xyz(line2)
                extend_vec = hydro_vec - ha_vec
                atom_hydro_vec[ha_label] = extend_vec
    #print(atom_hydro_vec)
    #print(atom_hydro_vec["C10"])

    candidate = sorted(set(candidate), key= candidate.index)
    #print(candidate)
    with open(ligand,'r')as lig1:
        for line in lig1:
            if (line[0:6]=="HETATM" or line[0:6]=="ATOM  "):
                if str(int(line[7:11])) in candidate:
                    with open(candidates,'a')as cdd:
                        print(line, end='', file=cdd)
    return candidates, atom_hydro_vec

