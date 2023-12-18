from script_sincho.basic.basic_func import *


def procheck(protein):
    ##################### make protein heavy atom dictionary #####################
    pro_lst = []
    for i in open(protein,'r').readlines():
        if i[0:6]=='ATOM  ':
            pro_lst.append(i)
    return pro_lst
    ##################### make protein heavy atom dictionary #####################
