from script_sincho.basic.basic_func import *
from script_sincho.extend_score.constant import *

def term3(pqr_inf, involve, set_atm, surround):

    ##################### constant value #####################
    atom_sym = fpocket_constant_value("atom_sym")
    electronegativity = fpocket_constant_value("electronegativity")
    vdw = fpocket_constant_value("vdw")
    ##################### constant value #####################

    ################ calculate X3 asa #######################
    asa_atm = 0
    tmp_num = 0
    for i in set_atm:
        #print(i)
        #revise 230607 kudo
        if i[76:78].strip().capitalize() in atom_sym:
            idx = atom_sym.index((i[76:78]).strip().capitalize())
            #idx = atom_sym.index(i.split()[-1].capitalize())

            atom = atom_sym[idx]
            elec = electronegativity[idx]
            r = vdw[idx]
        
            #set alpha selected#
            inv_alpha = []
            for k,v in involve.items():
                if i in v:
                    inv_alpha.append(k)

            if elec > 2.8:
                ax = float(i[30:38])
                ay = float(i[38:46])
                az = float(i[46:54])
                nop = 100
                inc = 3.1415926535897931*(3-math.sqrt(5))
                off = 2.0/nop
                nonburried = 0

                asa=0

                for k in range(nop):
                    py = k*off-1.0+(off/2.0)
                    rr = math.sqrt(1.0-py*py)
                    phi= k*inc
                    px = np.cos(phi)*rr
                    pz = np.sin(phi)*rr
                    x = ax+(r+2.2)*px
                    y = ay+(r+2.2)*py
                    z = az+(r+2.2)*pz

                    vec_p = np.array([x,y,z])
                    
                    # judge 1 #
                    vrefburried = 1
                    for a in inv_alpha:
                        if float(np.linalg.norm(vec_p - vec_xyz(a), ord=2))<=float(a[66:71]):
                            vrefburried = 0
                    # judge 2 #
                    burried = 0
                    for b in surround:
                        #revise 230607 kudo
                        if b[76:78].strip().capitalize() in atom_sym:
                            idx2 = atom_sym.index((b[76:78]).strip())
                            r2 = vdw[idx2]
                            #if i != b:
                            if float(np.linalg.norm(vec_p - vec_xyz(b), ord=2))<r2+2.2:
                                burried = 1
                                break
                    if vrefburried==0 and burried==0:
                        nonburried += 1
                        tmp_num += 1
                #calc. asa
                asa = (4*(3.1415926535897931)*(r+2.2)*(r+2.2)/nop)*(nonburried)
                #print(asa)
                asa_atm += asa
    ################ calculate X3 asa #######################

    return asa_atm
