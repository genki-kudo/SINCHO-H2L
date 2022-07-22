import numpy as np
from script.basic_func import lat_gen
from script.basic_func import t_file
import subprocess
import glob
import os
import shutil


def poc_vol_calc(outdirname):
    cls_list = glob.glob(outdirname+'/*.pqr')
    for cls in cls_list:
        path = os.path.splitext(cls)[0]
        name = (path.split('/')[-1])
        #print(name)
        lat_gen(cls, 'ATOM  ', outdirname+'/lat_'+name+'.pdb')

    estimate_pocvol = {}
    lat_list = glob.glob(outdirname+'/lat*')
    for lat in lat_list:
        with open(lat,'r')as file:
            estimate_pocvol[lat]=sum(1 for line in file)
    return estimate_pocvol



