from script_sincho.basic.basic_func import *

def involve_set(pqr, fpocket_output, pro_lst):
    ################# pqr_information and involved_atom information#####################
    pqr_inf = []#spheres_datas
    pqr_vec = []#spheres_vector
    pqr_rad = []#spheres_radius
    fpoc_original_poc_num = []#fpocket original pocket index
    involve = {}
    set_atm = []
    for i in open(pqr,'r').readlines():
        if i[0:6]=='ATOM  ' or i[0:6]=='HETATM':
            pqr_inf.append(i)
            pqr_vec.append(vec_xyz(i))
            pqr_rad.append(float(i[66:71]))
            #this check is dependent on fpocket's version. 
            # fpocket2->int(i[22:26])-1
            # fpocket4->int(i[22:26])
            if int(i[22:26]) not in fpoc_original_poc_num:
                fpoc_original_poc_num.append(int(i[22:26]))
    for i in range(len(pqr_inf)):
        inv = []
        for num in fpoc_original_poc_num:
            for j in open(fpocket_output+'/pockets/pocket'+str(num)+'_atm.pdb').readlines():
                if len(vec_xyz(j))==3:
                    if math.isclose(float(np.linalg.norm(pqr_vec[i] - vec_xyz(j))),float(pqr_rad[i]), rel_tol=0.005):
                        inv.append(j)
        involve[pqr_inf[i]]=inv
    #duplicate check
    for k,v in involve.items():
        for vi in v:
            if vi not in set_atm:
                set_atm.append(vi)
    surround = []
    for i in range(len(pqr_inf)):
        tmp_num = 0
        for j in pro_lst:
            if float(np.linalg.norm(pqr_vec[i]-vec_xyz(j),ord=2))<float(pqr_rad[i]+1.0):
                surround.append(j)
                tmp_num += 1
    ################# pqr_information and involved_atom information#####################

    return pqr_inf, involve, set_atm, surround
