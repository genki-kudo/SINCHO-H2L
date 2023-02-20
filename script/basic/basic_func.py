import subprocess
import math
import numpy as np
import pandas as pd
from Bio.PDB import *
from IPython.core.debugger import Pdb

def dict_pdb_noh(pdbfile):
    dict = {}
    for line in open(pdbfile).readlines():
        if (line[0:6]=="HETATM" or line[0:6]=="ATOM  ") and line[76:78]!=' H':
            dict[str(line[6:11])] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return dict

def dict_pdb_noh_nohetatm(pdbfile):
    dict = {}
    for line in open(pdbfile).readlines():
        if line[0:6]=="ATOM  " and line[76:78]!=' H':
            dict[str(line[6:11])] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return dict


#extract xyz of a atom
def ext_xyz(pdbline):
    if pdbline[0:6]=="HETATM" or pdbline[0:6]=="ATOM  ":
        x = float('{:.3f}'.format(float(pdbline[30:38])))
        y = float('{:.3f}'.format(float(pdbline[38:46])))
        z = float('{:.3f}'.format(float(pdbline[46:54])))
        return x, y, z
    else:
        return 'None'

#calculate vector of a atom
def vec_xyz(pdbline):
    if pdbline[0:6]=="HETATM" or pdbline[0:6]=="ATOM  ":
        x = float('{:.3f}'.format(float(pdbline[30:38])))
        y = float('{:.3f}'.format(float(pdbline[38:46])))
        z = float('{:.3f}'.format(float(pdbline[46:54])))
        xyz = [x, y, z]
        vec_xyz = np.array(xyz)
        return vec_xyz
    else:
        return 'None'
    
def dist_cf(vec_a, vec_b, d_min):
    distance = float(np.linalg.norm(vec_a - vec_b))
    if distance <= d_min:
        return distance
    else:
        return d_min

def dist_cf_max(vec_a, vec_b, d_max):
    distance = float(np.linalg.norm(vec_a - vec_b))
    if distance >= d_max:
        return distance
    else:
        return d_max

#truncated file  
def t_file(filename):
    with open(filename,'w')as file:
        file.truncate(0)


def appr(xxx):
    if xxx < 0:
        return '{:7.03f}'.format(xxx - 0.5)
    else:
        return '{:7.03f}'.format(xxx + 0.5)
        
def lat_gen(inputname, column, outputname):
    t_file('lat.pdb')
    p_num = 0
    with open(inputname,'r')as poc:
        for line in poc:
            if line[0:6] == column:
                p_num += 1
                num_pdb = '{:5}'.format(p_num)
                lxi = math.modf(float(line[30:38]))[1]
                lyi = math.modf(float(line[38:46]))[1]
                lzi = math.modf(float(line[46:54]))[1]
                llx = appr(lxi)
                lly = appr(lyi)
                llz = appr(lzi)
                with open('lat.pdb','a')as lat:
                    print('HETATM'+str(num_pdb)+'      PLA A   1      '+str(llx)+' '+str(lly)+' '+str(llz)+'  1.00 10.00           H', file=lat)
    t_file('lat_ex.txt')
    with open('lat.pdb','r')as poc:
        for line in poc:
            lat_x = float(line[30:38])
            lat_y = float(line[38:46])
            lat_z = float(line[46:54])
            for j in range(3):
                x = round(float(j-1.0),1)
                for k in range(3):
                    y = round(float(k-1.0),1)
                    for l in range(3):
                        z = round(float(l-1.0),1)
                        with open('lat_ex.txt','a')as exp:
                            print(lat_x+x, lat_y+y, lat_z+z, file=exp)
    df = pd.read_csv('lat_ex.txt', sep=" ",header=None)
    dup = df.drop_duplicates()
    t_file(outputname)
    p_num = 0
    for item in zip(dup[0],dup[1],dup[2]):
        p_num += 1
        num_pdb = '{:5}'.format(p_num)
        one_x = '{:7.03f}'.format(float(item[0]))
        one_y = '{:7.03f}'.format(float(item[1]))
        one_z = '{:7.03f}'.format(float(item[2]))
        with open (outputname,'a')as poc:
            print('HETATM'+num_pdb+'      PLA A   1     '+one_x+' '+one_y+' '+one_z+'  1.00 10.00           H', file=poc)
    subprocess.call(['rm', 'lat.pdb'])
    subprocess.call(['rm', 'lat_ex.txt'])

def Bio_position(file, num):
    parser = PDBParser()
    data = parser.get_structure("obj", file)
    for model in data.get_list():
        for chain in model.get_list():
            for residue in chain.get_list():
                residue.get_resname()
                atom = residue.get_list()[num]
                name = atom.get_name()
                vec = np.array(atom.get_vector())
    return name, vec

def getNearestValue(list, num):
    idx = np.abs(np.asarray(list)-num).argmin()
    return list[idx]

def listing_pdbcoordinate(pdbfile, listname):
    with open(pdbfile, 'r')as pdb:
        for line in pdb:
            if line[0:6]=='HETATM' or line[0:6]=='ATOM  ':
                listname.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return listname

def lattice_generate_from_list(inlistname, outlistname):
    outlistname = []
    for a in inlistname:
        x, y, z = float(appr(math.modf(a[0])[1])), float(appr(math.modf(a[1])[1])), float(appr(math.modf(a[2])[1]))
        for j in range(3):
            x1 = round(float(j-1.0),1)
            for k in range(3):
                y1 = round(float(k-1.0),1)
                for l in range(3):
                    z1 = round(float(l-1.0),1)
                    outlistname.append([x+x1, y+y1, z+z1])
    outlistname = [tuple(i) for i in outlistname]
    #print(outlistname, len(outlistname))
    outlistname = list(set(outlistname))
    outlistname = [list(i) for i in outlistname]
    #print(outlistname, len(outlistname))
    return outlistname

def coordinate_of_atom(pdbfile, atomname):
    with open(pdbfile, 'r')as inp:
        for line in inp:
            if line[12:16].replace(' ','') == atomname:
                vec = np.array([float(line[30:38]),float(line[38:46]), float(line[46:54])])
        
    return vec

def coordinate_bonding_hydrogen(pdbfile, heavyatom):
    bondh_coor = []
    with open(pdbfile, 'r')as inp:
        for line1 in inp:
            if line1[12:16].replace(' ','') == heavyatom:
                hanum = int(line1[6:11])
        cand = []
    with open(pdbfile, 'r')as inp:
        for line2 in inp:
            if line2[0:6]=='CONECT' and int(line2[6:11]) == hanum:
                l = line2.split()
                tmp = 0
                for a in l:
                    tmp += 1
                    if tmp>=3:
                        cand.append(a)
    for i in cand:
        with open(pdbfile, 'r')as inp:
            for line3 in inp:
                # if (line3[0:6]=='ATOM  ' or line3[0:6]=='HETATM') and int(line3[6:11])==int(i) and line3[12:16].replace(' ','')=='H':
                if (line3[0:6]=='ATOM  ' or line3[0:6]=='HETATM') and int(line3[6:11])==int(i) and line3.split()[-1]=='H':
                    bondh_coor.append(np.array([float(line3[30:38]),float(line3[38:46]), float(line3[46:54])]))
    return bondh_coor

def theta_calc(a,b,c):
    b_to_a = np.array(a) - np.array(b)
    b_to_c = np.array(c) - np.array(b)
    norm_b_to_a = np.linalg.norm(b_to_a)
    norm_b_to_c = np.linalg.norm(b_to_c)
    cos_abc = np.dot(b_to_a, b_to_c)/(norm_b_to_a*norm_b_to_c)
    theta = np.arccos(cos_abc)*180/np.pi
    return theta

def norm_calc(a,b):
    b_to_a = np.array(a) - np.array(b)
    norm_b_to_a = np.linalg.norm(b_to_a)
    return norm_b_to_a


def get_unique_list(seq):
    seen = []
    return [x for x in seq if x not in seen and not seen.append(x)]
    

def pdb_out(x,y,z,output_name, num):
    one_x = '{:7.03f}'.format(float(x))
    one_y = '{:7.03f}'.format(float(y))
    one_z = '{:7.03f}'.format(float(z))
    with open (output_name,'a')as poc:
        print('HETATM'+'{:5}'.format(num)+'      PLA A   1     '+one_x+' '+one_y+' '+one_z+'  1.00 10.00           H', file=poc)


def logo():
    print('        PPPPPPPPP/     22222222/      CCCCCCCC/')
    print('       PP/     PP/           22/     CC/       ')
    print('      PPPPPPPPP/     22222222/      CC/')
    print('     PP/            22/            CC/       ')
    print('    PP/            22222222/      CCCCCCCC/')
    print('convert sincho')




