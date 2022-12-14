from script.poc_vol_calc import poc_cbn_conv
from script.clus_mw_estimate import poc_mw_calc
from script.clus_mw_estimate import lig_poc_dist_list
from script.clus_mw_estimate import estimate_idealmw
from script.clus_mw_estimate import extend_startatom
from script.points_rest import r_candidate
import sys
import json
import os

#absolute path of script file 
absscriptpath = os.path.abspath(__file__)

#input
args = sys.argv
ligand = args[1]
ideal_mw = int(args[2])

#stat-correlation-set
mw_from_vol = json.loads(open(absscriptpath+'/statistics_dataset/vol_mw_fitting.json').read())
mw_from_dist = json.loads(open(absscriptpath+'/statistics_dataset/d_mw_fitting.json').read())

#output-setup
candidates = 'R-candidate.pdb'
clusdir = './p2c_output/cluster/'

#atom-selection's cosine
cos_value = 0

#estimate mw from pocket-volume
estimate_pocvol = poc_cbn_conv(clusdir)
estimate_poc_mw = poc_mw_calc(estimate_pocvol, mw_from_vol)

#select extend candidates based on the hydrogen and angle
ligand, atom_hydro_vec = r_candidate(ligand, candidates)

#calculate distance between atoms and pockets
dist_dict = lig_poc_dist_list(clusdir, ligand, atom_hydro_vec, cos_value)
#print(dist_dict)
#print(sorted(dist_dict.items(),key=lambda x:x[1]))
#print(len(dist_dict))

#estimate ideal mw based on pocket-volume and distance
extend_list = estimate_idealmw(estimate_poc_mw, ligand, dist_dict)

#print(extend_list)
print(sorted(extend_list.items(), key=lambda x:x[1]))
#print(extend_list["cluster0_C9 "])
for i in extend_list:
    if "cluster0_" in i:
        print(i)
        print(dist_dict[i])


goal1, goal2, goal3 = extend_startatom(extend_list, ideal_mw)
print("predicted next-pocket and start-atom: No.1->", goal1)
print("predicted next-pocket and start-atom: No.2->", goal2)
print("predicted next-pocket and start-atom: No.3->", goal3)
