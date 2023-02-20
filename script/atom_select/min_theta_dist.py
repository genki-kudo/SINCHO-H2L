from script.basic.basic_func import *

def min_theta_dist(candidate_atom_and_hydro, ligand, pqr):
    atom_theta_dist = []
    clusatom_coor = []
    listing_pdbcoordinate(pqr, clusatom_coor)
    for a, h in candidate_atom_and_hydro.items():
        theta_min = 180
        for bh in h:
            heavyatom_coor = coordinate_of_atom(ligand,a.split('_')[-1])
            theta_c = theta_calc(bh, heavyatom_coor, np.mean(clusatom_coor,axis=0))
            if theta_c <= theta_min:
                theta_min = theta_c
            dist_c = norm_calc(heavyatom_coor, np.mean(clusatom_coor,axis=0))
        atom_theta_dist.append([a, theta_min, dist_c])

    return atom_theta_dist

                