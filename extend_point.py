from script.volume_mw import vol_mw
from script.poc_vol_calc import poc_vol_calc
from script.clus_mw_estimate import poc_mw_calc
from script.clus_mw_estimate import lig_poc_dist_list
from script.clus_mw_estimate import estimate_idealmw
from script.clus_mw_estimate import extend_startatom
from script.points_rest import r_candidate
import sys

args = sys.argv

res1 = vol_mw('/home/user01/GIT-repository/02-EXTEND-POINT/dataset/')
ligand = args[1]
candidates = 'R-candidate.pdb'
clusdir = './p2c_output/cluster/'
cos_value = 0
ideal_mw = int(args[2])

estimate_pocvol = poc_vol_calc(clusdir)
estimate_poc_mw = poc_mw_calc(estimate_pocvol, res1)
#print(estimate_poc_mw)

ligand, atom_hydro_vec = r_candidate(ligand, candidates)

dist_dict = lig_poc_dist_list(clusdir, ligand, atom_hydro_vec, cos_value)
#print(dist_dict)
print(sorted(dist_dict.items(),key=lambda x:x[1]))
print(len(dist_dict))

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
